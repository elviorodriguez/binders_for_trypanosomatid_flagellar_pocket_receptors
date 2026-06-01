#!/usr/bin/env python
"""
This script takes two PDB files:
  - pdb_a: the physiological (original) model containing two chains:
       • a physiological target (e.g. chain A)
       • a physiological binder (e.g. chain B) with a conserved motif.
  - pdb_b: the RFdiffusion output model containing two chains:
       • a diffusion binder (e.g. chain A, all G's) 
       • a diffusion target (e.g. chain B, all G's)

It then:
  1. Applies a hard-coded rotation+translation (from a Matchmaker log) to both chains from pdb_a.
  2. Replaces the diffusion target with the transformed physiological target.
  3. In the diffusion binder, for each residue in the physiological binder that is part of the conserved motif 
     (specified by a list of residue numbers from pdb_a), finds the corresponding residue in the diffusion binder
     by comparing CA positions (if within a threshold) and grafts the side-chain atoms (and residue name)
     from the physiological binder onto the diffusion binder.
  4. Writes out a new PDB with the modified chains.
  
Usage:
  Single reference:
    python script.py pdb_a.pdb --pdb_b diff.pdb --phys_target A --phys_binder B --diff_binder A --diff_target B --motif_residues 28 30 31 --output_dir outdir

  Batch reference:
    python script.py pdb_a.pdb --batch_pdb_b ./diff_dir --phys_target A --phys_binder B --diff_binder A --diff_target B --motif_residues 28 30 31 --output_dir outdir

Notes:
  - The transformation (rotation+translation) is hard-coded.
  - The CA-based grafting uses a threshold (default 1.5 Å) to determine a match.
"""

import os
import argparse
from copy import deepcopy
import copy
import numpy as np
from Bio.PDB import PDBParser, PDBIO, is_aa

#############################
# Helper Functions
#############################

def load_structure(pdb_path, structure_id=None):
    """Load a structure from a PDB file."""
    parser = PDBParser(QUIET=True)
    if structure_id is None:
        structure_id = os.path.splitext(os.path.basename(pdb_path))[0]
    return parser.get_structure(structure_id, pdb_path)

def extract_chain(structure, chain_id):
    """Return a deep copy of the specified chain from the structure."""
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return chain.copy()
    raise ValueError(f"Chain {chain_id} not found in structure.")

def apply_transformation(chain, rotation_matrix, translation_vector):
    """Apply a rotation and translation to every atom in the chain (in-place)."""
    for atom in chain.get_atoms():
        coord = atom.get_coord()
        new_coord = np.dot(rotation_matrix, coord) + translation_vector
        atom.set_coord(new_coord)
    return chain

def graft_conserved_sidechains(phys_binder, diff_binder, conserved_resnums, threshold=1.5):
    """
    For each residue in the transformed physiological binder (phys_binder) whose residue number is in conserved_resnums,
    find the residue in the diffusion binder (diff_binder) that has the closest CA atom (if within threshold).
    Then, update the diffusion binder residue:
      - Set its residue name to that of the physiological binder residue.
      - Remove its side-chain atoms (all atoms except backbone: N, CA, C, O).
      - Copy over the side-chain atoms from the physiological binder residue (deep copy).
    
    Print debug information.
    """
    backbone = {"N", "CA", "C", "O"}
    # Build a mapping for diffusion binder: residue id -> CA coordinate
    diff_residues = []
    for res in diff_binder:
        if is_aa(res) and 'CA' in res:
            diff_residues.append(res)
    
    for res_phys in phys_binder:
        phys_resnum = res_phys.get_id()[1]
        if phys_resnum not in conserved_resnums:
            continue
        if 'CA' not in res_phys:
            continue
        phys_ca = res_phys['CA'].get_coord()
        # Find closest diffusion residue by CA distance:
        best_res = None
        best_dist = float('inf')
        for res_diff in diff_residues:
            if 'CA' not in res_diff:
                continue
            diff_ca = res_diff['CA'].get_coord()
            dist = np.linalg.norm(phys_ca - diff_ca)
            if dist < best_dist:
                best_dist = dist
                best_res = res_diff
        if best_res is None or best_dist > threshold:
            print(f"Warning: No matching residue found in diffusion binder for phys residue {phys_resnum} (best distance {best_dist:.2f}).")
            continue
        # Graft side-chain atoms: update residue name, remove non-backbone atoms, and copy side-chain atoms.
        old_name = best_res.resname
        best_res.resname = res_phys.resname
        print(f"Grafting: Phys residue {phys_resnum} ({res_phys.resname}) -> Diff residue {best_res.get_id()[1]} (was {old_name}), CA distance = {best_dist:.2f}")
        for atom in list(best_res):
            if atom.get_name() not in backbone:
                best_res.detach_child(atom.get_id())
        for atom in res_phys:
            if atom.get_name() not in backbone:
                new_atom = copy.deepcopy(atom)
                best_res.add(new_atom)
                print(f"  Added atom {new_atom.get_name()} to residue {best_res.get_id()[1]}.")

def replace_chain_in_structure(ref_structure, new_chain, chain_id):
    """
    In ref_structure (assumed to have one model), remove any chain with id chain_id and add new_chain (with chain_id set).
    Returns a deep copy of the modified structure.
    """
    structure_copy = deepcopy(ref_structure)
    model = structure_copy[0]
    if chain_id in model:
        model.detach_child(chain_id)
    new_chain.id = chain_id
    model.add(new_chain)
    return structure_copy

#############################
# Main Processing Function
#############################

def process_pair(pdb_a_path, pdb_b_path, phys_target_id, phys_binder_id,
                 diff_target_id, diff_binder_id, conserved_resnums,
                 rotation_matrix, translation_vector, output_dir, threshold):
    """
    Process one pair:
      - Load pdb_a (physiological model) and pdb_b (RFdiffusion output).
      - Extract from pdb_a the physiological target (phys_target_id) and physiological binder (phys_binder_id).
      - Extract from pdb_b the diffusion binder (diff_binder_id) and diffusion target (diff_target_id).
      - Apply the transformation (rotation+translation) to both chains from pdb_a.
      - Replace the diffusion target in pdb_b with the transformed physiological target.
      - In the diffusion binder, graft conserved side-chain atoms from the transformed physiological binder.
      - Write out the modified pdb_b as "mod_<basename>.pdb" in output_dir.
    """
    # Construct output file name.
    base_b = os.path.splitext(os.path.basename(pdb_b_path))[0]
    out_path = os.path.join(output_dir, f"mod_{base_b}.pdb")

    # Load structures.
    phys_struct = load_structure(pdb_a_path, "phys")
    diff_struct = load_structure(pdb_b_path, "diff")
    
    # Extract chains.
    phys_target = extract_chain(phys_struct, phys_target_id)
    phys_binder = extract_chain(phys_struct, phys_binder_id)
    diff_binder = extract_chain(diff_struct, diff_binder_id)
    diff_target = extract_chain(diff_struct, diff_target_id)
    
    # Apply the transformation (same for both chains).
    phys_target = apply_transformation(phys_target, rotation_matrix, translation_vector)
    phys_binder = apply_transformation(phys_binder, rotation_matrix, translation_vector)
    
    # Replace the diffusion target with the transformed physiological target.
    mod_struct = replace_chain_in_structure(diff_struct, phys_target, diff_target_id)
    
    # In the diffusion binder, graft conserved side chains from the physiological binder.
    graft_conserved_sidechains(phys_binder, diff_binder, conserved_resnums, threshold=threshold)
    # Replace the diffusion binder in the structure.
    mod_struct = replace_chain_in_structure(mod_struct, diff_binder, diff_binder_id)
    
    # Write out the modified structure.
    io = PDBIO()
    io.set_structure(mod_struct)
    io.save(out_path)
    print(f"Saved modified structure to {out_path}")

#############################
# Main CLI Entry Point
#############################

def main():
    parser = argparse.ArgumentParser(
        description="Graft conserved motif side-chain atoms from physiological binder to RFdiffusion binder, and replace target chain."
    )
    parser.add_argument("pdb_a", help="Path to physiological PDB file (with real side chains).")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pdb_b", help="Path to a single RFdiffusion output PDB file (with all G's).")
    group.add_argument("--batch_pdb_b", help="Directory with RFdiffusion output PDB files (first level only).")
    parser.add_argument("--phys_target", required=True, help="Chain ID in pdb_a for the physiological target.")
    parser.add_argument("--phys_binder", required=True, help="Chain ID in pdb_a for the physiological binder (with conserved motif).")
    parser.add_argument("--diff_target", required=True, help="Chain ID in pdb_b for the diffusion target (to be replaced).")
    parser.add_argument("--diff_binder", required=True, help="Chain ID in pdb_b for the diffusion binder (to receive grafted side chains).")
    parser.add_argument("--motif_residues", type=int, nargs='+', required=True,
                        help="List of conserved motif residue numbers from the physiological binder (pdb_a).")
    parser.add_argument("--output_dir", required=True, help="Directory to write modified structures.")
    parser.add_argument("--threshold", type=float, default=1.5, help="CA distance threshold (in Å) for matching residues (default: 1.5).")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Hard-coded transformation (from your Matchmaker log)
    rotation_matrix = np.array([
        [0.99999822, -0.00132335, 0.00134424],
        [0.00131753, 0.99998979, 0.00432183],
        [-0.00134995, -0.00432005, 0.99998976]
    ])
    translation_vector = np.array([9.84883114, -15.95180552, 5.08449227])
    
    if args.pdb_b:
        process_pair(
            pdb_a_path=args.pdb_a,
            pdb_b_path=args.pdb_b,
            phys_target_id=args.phys_target,
            phys_binder_id=args.phys_binder,
            diff_target_id=args.diff_target,
            diff_binder_id=args.diff_binder,
            conserved_resnums=args.motif_residues,
            rotation_matrix=rotation_matrix,
            translation_vector=translation_vector,
            output_dir=args.output_dir,
            threshold=args.threshold
        )
    else:
        for fname in os.listdir(args.batch_pdb_b):
            if fname.lower().endswith(".pdb"):
                pdb_b_path = os.path.join(args.batch_pdb_b, fname)
                process_pair(
                    pdb_a_path=args.pdb_a,
                    pdb_b_path=pdb_b_path,
                    phys_target_id=args.phys_target,
                    phys_binder_id=args.phys_binder,
                    diff_target_id=args.diff_target,
                    diff_binder_id=args.diff_binder,
                    conserved_resnums=args.motif_residues,
                    rotation_matrix=rotation_matrix,
                    translation_vector=translation_vector,
                    output_dir=args.output_dir,
                    threshold=args.threshold
                )

if __name__ == "__main__":
    main()

