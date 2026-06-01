#!/bin/bash

# This script contains the code that I have run to generate MRG binders using the scaffolding method

# Run RFdiffusion
/home/elvio/RFdiffusion/scripts/run_inference.py \
'inference.input_pdb=../input_models/MRG_MRGBP_truncated.pdb' \
'contigmap.contigs=[A425-461/A505-641/A652-683/0 5-50/B28-31/5-50]' \
'inference.num_designs=1000' \
'contigmap.inpaint_seq=[B29]' \
'contigmap.length=50-80' \
'ppi.hotspot_res=[A433,A436,A450,A446,A561,A564,A565,A566,A567,A568,A572,A623,A624,A626,A627,A630,A631,A634]' \
inference.ckpt_override_path=/home/elvio/RFdiffusion/models/Complex_base_ckpt.pt \
'potentials.guiding_potentials=["type:binder_ROG,weight:10,min_dist:5", "type:binder_ncontacts,weight:10", "type:interface_ncontacts,weight:5"]' \
potentials.guide_scale=2 \
potentials.guide_decay="quadratic"

# Transferring side-chains (here is needed to recover the translation vectore and rotation matrix using ChimeraX)
python 2-align_and_replace_with_motif_sidechain.py ../input_models/MRG_MRGBP_truncated_3.pdb --batch_pdb_b ../RFdiff_designs/samples/ --phys_target A --phys_binder B --diff_target A --diff_binder B --motif_residues 28 30 31 --output_dir ../pMPNN_designs/mpnn_input

# Run ProteinMPNN
./run_pMPNN_batch.sh -i ../pMPNN_designs/mpnn_input -o ../pMPNN_designs/mpnn_designs

# Generate the input for DiscobaMultimer
python 4-generate_discoba_multimer_input.py -i ../pMPNN_designs/mpnn_designs -t ../input_models/MRG.fasta -o ../discoba_multimer_designs -p ../input_models/MRG_truncated_isolated_target.pdb

# Run DiscobaMultimer
cd ../discoba_multimer_designs
discoba_multimer_batch -a -i merged_MSA -c AF2_templates.config database.fasta IDs_binders.txt | tee report.log

# Analyze desings
 python 6-analyze_binders.py --af_folder ../discoba_multimer_designs/AF2 --rf_folder ../RFdiff_designs/samples --output_dir ../analysis_plots

