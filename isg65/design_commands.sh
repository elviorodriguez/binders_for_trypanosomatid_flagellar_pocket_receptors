# cd to the directory to the project and work from there

# Set your own paths:
RFD_path=/home/esteban/RFdiffusion

# ---------------------- RFdiffusion ----------------------

# Activate RFdiff environment
conda activate SE3nv

# Run RFdiff for ISG65 binder based on C3dg motif (Motif scaffolding:sample_6)
$RFD_path/scripts/run_inference.py \
'inference.input_pdb=/home/esteban/Desktop/miniproteins_design/isg65/input/input_C3b_ISG65_clean.pdb' \
'contigmap.contigs=[A27-86/A94-154/A196-229/A251-316/0 5-30/B1164-1174/15-25/B1109-1114/5-30]' \
'inference.num_designs=10' \
'contigmap.inpaint_seq=[B1111,B1112,B1165,B1166,B1167,B1168,B1169,B1172,B1173]' \
'contigmap.length=50-80' \
'ppi.hotspot_res=[A70,A73,A74,A77,A282,A283,A286,A289,A290,A293]' \
inference.ckpt_override_path=$RFD_path/models/Complex_base_ckpt.pt \
'potentials.guiding_potentials=["type:binder_ROG,weight:5,min_dist:5", "type:binder_ncontacts,weight:10", "type:interface_ncontacts,weight:5"]' \
potentials.guide_scale=0.5 \
potentials.guide_decay="linear"

# ---------------------- ProteinMPNN ----------------------

# Run ProteinMPNN (ProteinMPNN environment is activated within run_pMPNN_batch.sh)
../scripts/3-run_pMPNN_batch.sh -i samples -o output_pMPNN

# ---------------------- AF2 ----------------------

# Desactivate any env and activate env to run AF2-multimer (DiscobaMultimer)
conda deactivate; conda activate DiscobaMultimer

# Generate the input for DiscobaMultimer (implementation of LocalColabFold)
python ../scripts/4-generate_discoba_multimer_input.py -i output_pMPNN -t input/isg65_7IP6.fa -o discoba_multimer_designs -p input/input_C3b_ISG65_clean.pdb

# The latter command generated a directory. cd to it
cd discoba_multimer_designs

# Execute AF2
discoba_multimer_batch -a -i merged_MSA -c AF2_templates.config database.fasta IDs_binders.txt | tee report.log

# ---------------------- Score models ----------------------

# Deactivate environment and activate one that can permita manipulate tables and PDB structures with BioPython (here we use MultimerMapper)
conda deactivate; conda activate MultimerMapper

# cd to the project directory once more

# Analyze desings
python ../scripts/6-analyze_binders.py --af_folder discoba_multimer_designs/AF2 --rf_folder samples --output_dir analysis_plots
