#!/usr/bin/env python3
"""
Compile designed protein binder sequences into a unified database
and create pairwise A3M alignments of binders with one or more target chains.

This script processes all FASTA files in the ProteinMPNN output directories and
compiles them into a single database.fasta file with standardized sequence headers.
It also creates:
  - IDs_binders.txt  : TSV mapping every target chain ID + binder ID (one line per binder)
  - merged_MSA/*.a3m : Pairwise A3M alignments suitable for ColabFold / AF2 multimer input

The target FASTA may contain ONE sequence (monomer) or MULTIPLE sequences (oligomer).
All sequences found in the target FASTA are treated as fixed target chains.

Usage:
    python compile_binders.py -i INPUT_DIR -t TARGET_FASTA -p PDB_TEMPLATE [-o OUTPUT_DIR]

Required arguments:
    -i, --input_dir     Directory containing ProteinMPNN output folders (from design_binders.sh)
    -t, --target_fasta  FASTA file with one or more target chain sequences
    -p, --pdb_template  PDB file used as structural template (binder + target complex)

Optional arguments:
    -o, --output_dir    Directory to save output files (default: current directory)
    -h, --help          Show this help message and exit

Output files:
    database.fasta        : All target sequences + all designed binder sequences
    IDs_binders.txt       : <target1_ID> [<target2_ID> ...] <binder_ID>  (one line per binder)
    merged_MSA/*.a3m      : Per-binder A3M files with paired + unpaired rows for every chain
    templates/targ.pdb    : Copy of the template PDB
    AF2_templates.config  : Ready-to-use AF2 config file
"""

import os
import sys
import glob
import argparse
import shutil
import re


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Compile ProteinMPNN designed sequences into a unified database "
            "and create pairwise A3M alignments for AF2 multimer prediction."
        )
    )
    parser.add_argument("-i", "--input_dir",    required=True,
                        help="Directory containing ProteinMPNN output folders")
    parser.add_argument("-t", "--target_fasta", required=True,
                        help="FASTA file with one or more target chain sequences")
    parser.add_argument("-o", "--output_dir",   default=".",
                        help="Directory to save output files (default: current directory)")
    parser.add_argument("-p", "--pdb_template", required=True,
                        help="PDB file used as structural template (binder + target complex)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------

def read_target_fasta(target_fasta_path):
    """
    Read one or more sequences from a FASTA file.

    Returns
    -------
    list of (header: str, sequence: str)
        One entry per '>' record found in the file, in order.
    """
    if not os.path.exists(target_fasta_path):
        raise FileNotFoundError(f"Target FASTA file not found: {target_fasta_path}")

    records = []
    current_header = None
    current_seq_parts = []

    with open(target_fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save previous record before starting a new one
                if current_header is not None:
                    records.append((current_header, "".join(current_seq_parts)))
                current_header = line.lstrip(">")
                current_seq_parts = []
            else:
                current_seq_parts.append(line)

    # Save last record
    if current_header is not None:
        records.append((current_header, "".join(current_seq_parts)))

    if not records:
        raise ValueError(f"No FASTA records found in: {target_fasta_path}")

    return records   # [(header, sequence), ...]


# ---------------------------------------------------------------------------
# ProteinMPNN output parsing
# ---------------------------------------------------------------------------

def process_binder_fastas(input_dir):
    """
    Collect all designed sequences from ProteinMPNN output directories.

    Returns
    -------
    list of (design_id: int, sample_id: int, sequence: str)
        Only the designed (sampled) sequences — the backbone/template sequence
        that ProteinMPNN writes as the first record is skipped.
    """
    binder_sequences = []
    design_dirs = sorted(
        d for d in os.listdir(input_dir)
        if os.path.isdir(os.path.join(input_dir, d))
    )

    for design_dir in design_dirs:
        design_id = int(design_dir.split("_")[-1])
        seq_dir = os.path.join(input_dir, design_dir, "seqs")

        if not os.path.exists(seq_dir):
            print(f"Warning: No 'seqs' directory found in {os.path.join(input_dir, design_dir)}")
            continue

        fasta_files = (
            glob.glob(os.path.join(seq_dir, "*.fa")) +
            glob.glob(os.path.join(seq_dir, "*.fasta"))
        )
        if not fasta_files:
            print(f"Warning: No FASTA files found in {seq_dir}")
            continue

        fasta_path = fasta_files[0]

        with open(fasta_path) as fh:
            lines = fh.readlines()

        current_seq = ""
        sample_id = 0

        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                if "sample=" in line:
                    if current_seq and sample_id > 0:
                        binder_sequences.append((design_id, sample_id, current_seq))
                    match = re.search(r"sample=(\d+)", line)
                    sample_id = int(match.group(1)) if match else sample_id + 1
                    current_seq = ""
                else:
                    # Backbone template record — skip
                    current_seq = ""
                    sample_id = 0
            elif sample_id > 0:
                current_seq += line

        if current_seq and sample_id > 0:
            binder_sequences.append((design_id, sample_id, current_seq))

    return binder_sequences


# ---------------------------------------------------------------------------
# Output file writers
# ---------------------------------------------------------------------------

def write_database_files(target_records, binder_sequences, output_dir):
    """
    Write database.fasta and IDs_binders.txt.

    database.fasta layout
    ---------------------
    >target1_header
    <sequence>
    [>target2_header   (if dimer / oligomer)
    <sequence>]
    ...
    >design_<N>.<M>
    <sequence>
    ...

    IDs_binders.txt layout (tab-separated, one line per binder)
    ---------------------
    target1_ID  [target2_ID  ...]  binder_ID

    Parameters
    ----------
    target_records : list of (header, sequence)
    binder_sequences : list of (design_id, sample_id, sequence)
    output_dir : str
    """
    os.makedirs(output_dir, exist_ok=True)

    fasta_path = os.path.join(output_dir, "database.fasta")
    tsv_path   = os.path.join(output_dir, "IDs_binders.txt")

    # --- database.fasta ---
    with open(fasta_path, "w") as fh:
        # All target chains first
        for header, seq in target_records:
            fh.write(f">{header}\n{seq}\n")

        # All designed binders
        for design_id, sample_id, seq in binder_sequences:
            binder_id = f"design_{design_id}.{sample_id}"
            fh.write(f">{binder_id}\n{seq}\n")

    # --- IDs_binders.txt ---
    target_ids = [header for header, _ in target_records]

    with open(tsv_path, "w") as fh:
        for design_id, sample_id, _ in binder_sequences:
            binder_id = f"design_{design_id}.{sample_id}"
            # tab-separated: target1_ID  [target2_ID ...]  binder_ID
            fields = [binder_id] + target_ids
            fh.write("\t".join(fields) + "\n")

    return fasta_path, tsv_path


def create_pairwise_a3m_alignments(binder_sequences, target_records, output_dir):
    """
    Create one A3M file per designed binder, suitable for AF2 / ColabFold multimer.

    A3M file layout (example: 2 target chains + 1 binder = 3 chains total)
    -----------------------------------------------------------------------
    #{binder_len},{target1_len},{target2_len}\\t1,1,1
    >101\\t102\\t103
    {binder_seq}{target1_seq}{target2_seq}          <- fully paired row
    >101
    {binder_seq}{'-'*target1_len}{'-'*target2_len}  <- binder unpaired
    >102
    {'-'*binder_len}{target1_seq}{'-'*target2_len}  <- target1 unpaired
    >103
    {'-'*(binder_len+target1_len)}{target2_seq}      <- target2 unpaired

    The chain ordering is: binder first, then targets in the order they appear
    in the target FASTA. This matches the chain order expected by AF2 when the
    binder is chain C and the targets are chains A and B (AF2 processes chains
    in the order they are presented in the MSA / sequence input).

    Parameters
    ----------
    binder_sequences : list of (design_id, sample_id, sequence)
    target_records   : list of (header, sequence)
    output_dir       : str

    Returns
    -------
    str : path to the merged_MSA directory
    """
    merged_msa_dir = os.path.join(output_dir, "merged_MSA")
    os.makedirs(merged_msa_dir, exist_ok=True)

    target_headers = [h for h, _ in target_records]
    target_seqs    = [s for _, s in target_records]
    target_lens    = [len(s) for s in target_seqs]
    n_targets      = len(target_records)

    # Build a short label for the filename that encodes all target IDs
    targets_label = "__vs__".join(target_headers)

    for design_id, sample_id, binder_seq in binder_sequences:
        binder_id  = f"design_{design_id}.{sample_id}"
        binder_len = len(binder_seq)

        # All chain lengths in order: binder, target1, target2, ...
        all_lens = [binder_len] + target_lens

        # --- Header line ---
        # Format: #{len1},{len2},...\t{stoich1},{stoich2},...
        lengths_str    = ",".join(str(l) for l in all_lens)
        stoichiometry  = ",".join(["1"] * (1 + n_targets))
        header_line    = f"#{lengths_str}\t{stoichiometry}"

        # --- Chain IDs for the A3M records ---
        # Convention: 101 = binder, 102 = target1, 103 = target2, ...
        binder_chain_id  = "101"
        target_chain_ids = [str(102 + i) for i in range(n_targets)]

        # Paired row: all chains concatenated
        paired_seq    = binder_seq + "".join(target_seqs)
        paired_id_str = "\t".join([binder_chain_id] + target_chain_ids)

        # Unpaired rows: one per chain, dashes everywhere else
        #   binder unpaired  → binder_seq + dashes for every target
        #   target_i unpaired → dashes for binder + dashes for targets before i
        #                        + target_i_seq + dashes for targets after i

        # Compute cumulative lengths for positioning dashes
        # Position of each chain in the concatenated sequence:
        #   binder  : [0 : binder_len]
        #   target0 : [binder_len : binder_len + target_lens[0]]
        #   target1 : [binder_len + target_lens[0] : binder_len + target_lens[0] + target_lens[1]]
        #   etc.
        total_len = sum(all_lens)

        # Binder unpaired: binder_seq then dashes for all targets
        binder_unpaired = binder_seq + "-" * sum(target_lens)

        # Target unpaired rows
        target_unpaired_seqs = []
        for i, (t_seq, t_len) in enumerate(zip(target_seqs, target_lens)):
            prefix_dash_len = binder_len + sum(target_lens[:i])
            suffix_dash_len = total_len - prefix_dash_len - t_len
            target_unpaired_seqs.append("-" * prefix_dash_len + t_seq + "-" * suffix_dash_len)

        # --- Write A3M file ---
        a3m_filename = f"{binder_id}__vs__{targets_label}.a3m"
        a3m_path     = os.path.join(merged_msa_dir, a3m_filename)

        with open(a3m_path, "w") as fh:
            fh.write(header_line + "\n")

            # Fully paired row (binder + all targets)
            fh.write(f">{paired_id_str}\n{paired_seq}\n")

            # Binder unpaired
            fh.write(f">{binder_chain_id}\n{binder_unpaired}\n")

            # Each target unpaired
            for chain_id, unpaired_seq in zip(target_chain_ids, target_unpaired_seqs):
                fh.write(f">{chain_id}\n{unpaired_seq}\n")

    return merged_msa_dir


def create_templates_dir(output_dir, pdb_file):
    """
    Copy the template PDB into a templates/ subdirectory and write an AF2 config.

    Returns
    -------
    str : path to the templates directory
    """
    templates_dir = os.path.join(output_dir, "templates")
    os.makedirs(templates_dir, exist_ok=True)

    copied_pdb_path = shutil.copy(pdb_file, templates_dir)
    renamed_pdb_path = os.path.join(templates_dir, "targ.pdb")
    os.rename(copied_pdb_path, renamed_pdb_path)

    config_path = os.path.join(output_dir, "AF2_templates.config")
    config_contents = (
        "--num-models 5\n"
        "--num-recycle 3\n"
        "--rank iptm\n"
        "--templates\n"
        "--custom-template-path ./templates\n"
    )
    with open(config_path, "w") as cfg:
        cfg.write(config_contents)

    return templates_dir


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_arguments()

    try:
        # --- Read target sequences (one or more chains) ---
        print(f"Reading target sequence(s) from {args.target_fasta}...")
        target_records = read_target_fasta(args.target_fasta)
        n_targets = len(target_records)
        print(f"Found {n_targets} target chain(s): {[h for h, _ in target_records]}")

        # --- Collect designed binder sequences ---
        print(f"Processing binder sequences from {args.input_dir}...")
        binder_sequences = process_binder_fastas(args.input_dir)
        print(f"Found {len(binder_sequences)} designed binder sequences")

        # --- Write database.fasta and IDs_binders.txt ---
        print(f"Writing output files to {args.output_dir}...")
        fasta_path, tsv_path = write_database_files(
            target_records, binder_sequences, args.output_dir
        )

        # --- Create A3M alignments ---
        print("Creating pairwise A3M alignments...")
        merged_msa_dir = create_pairwise_a3m_alignments(
            binder_sequences, target_records, args.output_dir
        )

        # --- Copy template PDB and write AF2 config ---
        templates_dir = create_templates_dir(args.output_dir, args.pdb_template)

        # --- Summary ---
        print(f"\nDone.")
        print(f"  database.fasta       : {fasta_path}")
        print(f"  IDs_binders.txt      : {tsv_path}")
        print(f"  merged_MSA/          : {merged_msa_dir}  ({len(binder_sequences)} A3M files)")
        print(f"  templates/           : {templates_dir}")
        print(f"  AF2_templates.config : {os.path.join(args.output_dir, 'AF2_templates.config')}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
