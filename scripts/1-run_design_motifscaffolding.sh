#!/bin/bash

# Using full length MRG
/home/elvio/RFdiffusion/scripts/run_inference.py 'inference.input_pdb=../input_models/MRG_MRGBP.pdb' 'contigmap.contigs=[A425-699/0 10-40/B28-31/10-40]' 'inference.num_designs=10' 'contigmap.inpaint_seq=[B29]'

# Using truncated MRG and hotspots
/home/elvio/RFdiffusion/scripts/run_inference.py 'inference.input_pdb=../input_models/MRG_MRGBP_truncated.pdb' 'contigmap.contigs=[A425-461/A505-641/A652-683/0 10-40/B28-31/10-40]' 'inference.num_designs=10' 'contigmap.inpaint_seq=[B29]' 'ppi.hotspot_res=[A443,A446,A447,A565,A566,A568,A627,A630,A634]'

# Using truncated MRG, hotspots and modified potentials
/home/elvio/RFdiffusion/scripts/run_inference.py 'inference.input_pdb=../input_models/MRG_MRGBP_truncated.pdb' 'contigmap.contigs=[A425-461/A505-641/A652-683/0 20-40/B28-31/20-40]' 'inference.num_designs=10' 'contigmap.inpaint_seq=[B29]' 'ppi.hotspot_res=[A443,A446,A447,A565,A566,A568,A627,A630,A634]' 'potentials.guiding_potentials=["type:monomer_ROG,weight:1,min_dist:5"]' potentials.guide_scale=2 potentials.guide_decay="quadratic" inference.ckpt_override_path=/home/elvio/RFdiffusion/models/Complex_base_ckpt.pt


# This works really well
/home/elvio/RFdiffusion/scripts/run_inference.py \
'inference.input_pdb=../input_models/MRG_MRGBP_truncated.pdb' \
'contigmap.contigs=[A425-461/A505-641/A652-683/0 5-40/B28-31/5-40]' \
'inference.num_designs=5' \
'contigmap.inpaint_seq=[B29]' \
'contigmap.length=50-80' \
'ppi.hotspot_res=[A433,A436,A438,A443,A450,A446,A447,A561,A564,A565,A566,A567,A568,A572,A623,A624,A626,A627,A630,A631,A634]' \
inference.ckpt_override_path=/home/elvio/RFdiffusion/models/Complex_base_ckpt.pt \
'potentials.guiding_potentials=["type:binder_ROG,weight:10,min_dist:5", "type:binder_ncontacts,weight:10", "type:interface_ncontacts,weight:5"]' \
potentials.guide_scale=2 \
potentials.guide_decay="quadratic"

# Reducing the hotspot area (1000 models)
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

