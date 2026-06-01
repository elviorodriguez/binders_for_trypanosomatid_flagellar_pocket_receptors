# cd to the directory to the project and work from there

# Set your own paths to the different softwares:
RFD_path=/home/esteban/RFdiffusion

# ---------------------- RFdiffusion ----------------------

# Activate RFdiff environment
conda activate SE3nv

# Run RFdiff for ISG65 binder based on C3dg motif (Motif scaffolding:samples)
$RFD_path/scripts/run_inference.py \
'inference.input_pdb=/home/esteban/Desktop/miniproteins_design/hphbr/input/hphbr_clean.pdb' \
'contigmap.contigs=[D38-297/0 10-20/B92-92/3-3/B96-96/13-30/B41-45/10-20]' \
'inference.num_designs=1000' \
'contigmap.inpaint_seq=[B43]' \
'contigmap.length=50-80' \ 'ppi.hotspot_res=[D157,D160,D161,D164,D200,D201,D203,D56,D57,D59,D60,D63,D66,D67,D70,D64,D61,D204,D153,D68,D71]' \
inference.ckpt_override_path=$RFD_path/models/Complex_base_ckpt.pt \
'potentials.guiding_potentials=["type:binder_ROG,weight:5,min_dist:5", "type:binder_ncontacts,weight:20", "type:interface_ncontacts,weight:5"]' \
potentials.guide_scale=0.5 \
potentials.guide_decay="linear"

# ---------------------- ProteinMPNN ----------------------

# Run ProteinMPNN (ProteinMPNN environment is activated within run_pMPNN_batch.sh)
../scripts/3-run_pMPNN_batch.sh -i samples -o output_pMPNN -t "D"

# ---------------------- AF2 ----------------------

# Desactivate any env and activate env to run AF2-multimer (DiscobaMultimer)
conda deactivate; conda activate DiscobaMultimer

# Generate the input for DiscobaMultimer (implementation of LocalColabFold)
python ../scripts/4-generate_discoba_multimer_input.py -i output_pMPNN -t input/hphbr.fasta -o discoba_multimer_designs -p input/hphbr_clean.pdb

# The latter command generated a directory. cd to it
cd discoba_multimer_designs

# Execute AF2
discoba_multimer_batch -a -i merged_MSA -c AF2_templates.config database.fasta IDs_binders.txt | tee report.log

# ---------------------- Score models ----------------------


# Deactivate environment and activate one that can permita manipulate tables and strures (here we use MultimerMapper)
conda deactivate; conda activate MultimerMapper

# cd to the project directory once more

# Analyze desings
python ../scripts/6b-analyze_binders.py --af_folder discoba_multimer_designs/AF2 --rf_folder samples --output_dir analysis_plots --binder_chain_rf "B" --binder_chain_af "A"
