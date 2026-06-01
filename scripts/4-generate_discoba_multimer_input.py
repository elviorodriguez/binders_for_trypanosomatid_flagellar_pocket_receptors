#!/usr/bin/env python3
"""
Compile designed protein binder sequences into a unified database
and create pairwise A3M alignments of binders with target

This script processes all fasta files in the ProteinMPNN output directories and compiles them
into a single database.fasta file with standardized sequence headers. It also creates a TSV
file mapping the target sequence to each designed binder sequence. Additionally, it creates 
pairwise A3M alignments of each binder with the target sequence.

Usage:
    python script.py [-h] -i INPUT_DIR -t TARGET_FASTA [-o OUTPUT_DIR] -p PDB_TEMPLATE

Required arguments:
    -i, --input_dir     Directory containing ProteinMPNN output folders (from design_binders.sh)
    -t, --target_fasta  FASTA file containing the target protein sequence
    -p, --pdb_template  PDB file used as template (Initial binder+target PDB)

Optional arguments:
    -o, --output_dir    Directory to save output files (default: current directory)
    -h, --help          Show this help message and exit

Output files:
    - database.fasta: Combined FASTA file with the target sequence and all designed binder sequences
    - IDs_binders.txt: TSV file mapping target sequence ID to each binder sequence ID
    - merged_MSA/*.a3m: Pairwise A3M alignments of each binder with the target sequence
    - templates/: Directory containing the template PDB file for AF2
    - AF2_templates.config: Configuration file for running AF2 with templates
"""

import os
import sys
import glob
import argparse
import shutil
import re


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Compile ProteinMPNN designed sequences into a unified database and create pairwise A3M alignments"
    )
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="Directory containing ProteinMPNN output folders (from run_pMPNN_batch.sh)"
    )
    parser.add_argument(
        "-t", "--target_fasta",
        required=True,
        help="FASTA file containing the target protein sequence"
    )
    parser.add_argument(
        "-o", "--output_dir",
        default=".",
        help="Directory to save output files (default: current directory)"
    )
    parser.add_argument(
        "-p", "--pdb_template",
        required=True,
        help="PDB file used as template (Initial binder+target PDB)."
    )
    return parser.parse_args()


def read_target_fasta(target_fasta_path):
    """
    Read the target protein sequence from a FASTA file.
    
    Args:
        target_fasta_path (str): Path to the target FASTA file
        
    Returns:
        tuple: (header, sequence) of the target protein
    """
    # Check if file exists
    if not os.path.exists(target_fasta_path):
        raise FileNotFoundError(f"Target FASTA file not found: {target_fasta_path}")
    
    # Read the FASTA file
    with open(target_fasta_path, 'r') as f:
        lines = f.readlines()
    
    # Extract header and sequence
    if not lines or not lines[0].startswith('>'):
        raise ValueError(f"Invalid FASTA format in file: {target_fasta_path}")
    
    header = lines[0].strip().lstrip('>')
    sequence = ''.join(line.strip() for line in lines[1:])
    
    return header, sequence


def process_binder_fastas(input_dir):
    """
    Process all binder FASTA files in the ProteinMPNN output directories.
    
    Args:
        input_dir (str): Directory containing ProteinMPNN output folders
        
    Returns:
        list: List of tuples (design_id, sample_id, sequence)
    """
    binder_sequences = []
    design_dirs = sorted([d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))])
    
    for design_dir in design_dirs:

        # Extract design ID
        design_id = int(design_dir.split("_")[-1])
        
        # Look for the fasta file in the seqs subdirectory
        design_path = os.path.join(input_dir, design_dir)
        seq_dir = os.path.join(design_path, "seqs")
        
        if not os.path.exists(seq_dir):
            print(f"Warning: No 'seqs' directory found in {design_path}")
            continue
        
        # Find the FASTA file
        fasta_files = glob.glob(os.path.join(seq_dir, "*.fa")) + glob.glob(os.path.join(seq_dir, "*.fasta"))
        if not fasta_files:
            print(f"Warning: No FASTA files found in {seq_dir}")
            continue
        
        fasta_path = fasta_files[0]  # Take the first fasta file if multiple exist
        
        # Parse the FASTA file
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
        
        # Extract sequences (skip the first sequence which is the backbone)
        current_seq = ""
        sample_id = 0
        
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                # If this is a designed sequence (contains 'sample=' in the header)
                if 'sample=' in line:
                    # If we've already collected a sequence, save it
                    if current_seq and sample_id > 0:
                        binder_sequences.append((design_id, sample_id, current_seq))
                    # Extract the sample number
                    match = re.search(r'sample=(\d+)', line)
                    if match:
                        sample_id = int(match.group(1))
                    else:
                        sample_id += 1
                    current_seq = ""
                else:
                    # This is the backbone sequence, reset
                    current_seq = ""
                    sample_id = 0
            elif sample_id > 0:
                # Collect the sequence
                current_seq += line
        
        # Don't forget the last sequence
        if current_seq and sample_id > 0:
            binder_sequences.append((design_id, sample_id, current_seq))
    
    return binder_sequences


def write_database_files(target_header, target_sequence, binder_sequences, output_dir):
    """
    Write the compiled database.fasta and IDs_binders.txt files.
    
    Args:
        target_header (str): Header of the target sequence
        target_sequence (str): Sequence of the target protein
        binder_sequences (list): List of tuples (design_id, sample_id, sequence)
        output_dir (str): Directory to save output files
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Paths for output files
    fasta_path = os.path.join(output_dir, "database.fasta")
    tsv_path = os.path.join(output_dir, "IDs_binders.txt")
    
    # Write database.fasta
    with open(fasta_path, 'w') as f:
        # Write target sequence first
        f.write(f">{target_header}\n{target_sequence}\n")
        
        # Write all binder sequences
        for design_id, sample_id, sequence in binder_sequences:
            binder_id = f"design_{design_id}.{sample_id}"
            f.write(f">{binder_id}\n{sequence}\n")
    
    # Write IDs_binders.txt
    with open(tsv_path, 'w') as f:
        for design_id, sample_id, _ in binder_sequences:
            binder_id = f"design_{design_id}.{sample_id}"
            f.write(f"{binder_id}\t{target_header}\n")
    
    return fasta_path, tsv_path


def create_pairwise_a3m_alignments(binder_sequences, target_header, target_sequence, output_dir):
    """
    Create pairwise A3M alignments of each binder with the target sequence.
    
    Args:
        binder_sequences (list): List of tuples (design_id, sample_id, sequence)
        target_header (str): Header of the target protein
        target_sequence (str): Sequence of the target protein
        output_dir (str): Directory to save output files
        
    Returns:
        str: Path to the directory containing the A3M alignments
    """
    # Create merged_MSA directory
    merged_msa_dir = os.path.join(output_dir, "merged_MSA")
    os.makedirs(merged_msa_dir, exist_ok=True)
    
    # Process each binder sequence
    for design_id, sample_id, binder_seq in binder_sequences:
        binder_id = f"design_{design_id}.{sample_id}"
        
        # Calculate lengths
        binder_length = len(binder_seq)
        target_length = len(target_sequence)
        
        # Create output A3M file path
        output_a3m_file = os.path.join(merged_msa_dir, f"{binder_id}__vs__{target_header}.a3m")
        
        with open(output_a3m_file, 'w') as f:
            # Write header line with dimensions
            f.write(f"#{binder_length},{target_length}\t1,1\n")
            
            # Write first sequence (concatenated binder + target)
            f.write(f">101\t102\n{binder_seq}{target_sequence}\n")
            
            # Write binder sequence with dashes for target
            f.write(f">101\n{binder_seq}{'-' * target_length}\n")
            
            # Write target sequence with dashes for binder
            f.write(f">102\n{'-' * binder_length}{target_sequence}")
            
            # Note: No longer include sequences from an input A3M file
    
    return merged_msa_dir


def create_templates_dir(output_dir, pdb_file):
    """
    Create templates directory and copy/rename the template PDB file.
    
    Args:
        output_dir (str): Directory to save output files
        pdb_file (str): Path to the template PDB file
        
    Returns:
        str: Path to the templates directory
    """
    # Create templates directory
    templates_dir = os.path.join(output_dir, "templates")
    os.makedirs(templates_dir, exist_ok=True)

    # Copy pdb file into templates directory and rename it
    copied_pdb_path = shutil.copy(pdb_file, templates_dir)    
    renamed_pdb_path = os.path.join(templates_dir, "targ.pdb")
    os.rename(copied_pdb_path, renamed_pdb_path)

    # Write AF2_templates.config
    config_path = os.path.join(output_dir, "AF2_templates.config")
    config_contents = '--num-models 5\n--num-recycle 3\n--rank iptm\n--templates\n--custom-template-path ./templates\n'
    with open(config_path, 'w') as cfg:
        cfg.write(config_contents)

    return templates_dir


def main():
    """Main function to process and compile binder sequences."""
    args = parse_arguments()
    
    try:
        # Read target sequence
        print(f"Reading target sequence from {args.target_fasta}...")
        target_header, target_sequence = read_target_fasta(args.target_fasta)
        
        # Process all binder sequences
        print(f"Processing binder sequences from {args.input_dir}...")
        binder_sequences = process_binder_fastas(args.input_dir)
        print(f"Found {len(binder_sequences)} designed binder sequences across multiple designs")
        
        # Write output files
        print(f"Writing output files to {args.output_dir}...")
        fasta_path, tsv_path = write_database_files(
            target_header, target_sequence, binder_sequences, args.output_dir
        )
        
        # Create pairwise A3M alignments
        print(f"Creating pairwise A3M alignments...")
        merged_msa_dir = create_pairwise_a3m_alignments(
            binder_sequences, target_header, target_sequence, args.output_dir
        )
        
        # Create templates dir and cfg
        templates_dir = create_templates_dir(args.output_dir, args.pdb_template)
        
        print(f"Successfully created database FASTA: {fasta_path}")
        print(f"Successfully created binder mapping TSV: {tsv_path}")
        print(f"Successfully created pairwise A3M alignments in: {merged_msa_dir}")
        print(f"Successfully created templates directory in: {templates_dir}")
        print(f"Created {len(binder_sequences)} A3M alignment files")
        
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
