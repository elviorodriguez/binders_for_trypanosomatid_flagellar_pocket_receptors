#!/bin/bash

conda activate MultimerMapper

# Without transferring side-chains
python align_and_replace.py MRG_MRGBP_truncated_2.pdb --batch_pdb_b ./samples_test --chain_a A --chain_b B --output_dir outdir

# Transferring side-chains
python align_and_replace_with_motif_sidechain.py MRG_MRGBP_truncated_2.pdb --batch_pdb_b samples_test/ --phys_target A --phys_binder B --diff_target B --diff_binder A --motif_residues 28 30 31 --output_dir outdir

