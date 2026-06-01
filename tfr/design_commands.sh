# cd to the directory to the project and work from there

# Set your own paths to the different softwares:
RFD_path=/home/esteban/RFdiffusion

# ---------------------- RFdiffusion ----------------------

# Activate RFdiff environment
conda activate SE3nv

# Run RFdiff for TbTfR (ESAG6/7) binder based on Tf motif (Motif scaffolding:sample_???)
$RFD_path/scripts/run_inference.py \
'inference.input_pdb=/home/esteban/Desktop/miniproteins_design/tfr/input/input_structure_clean.pdb' \
'contigmap.contigs=[A20-342/0 5-30/C349-370/5-30 B18-337/0]' \
'inference.num_designs=10' \
'contigmap.inpaint_seq=[C350,C351,C353,C354,C357,C358,C361,C362,C363,C365,C366,C368,C369,C370]' \
'contigmap.length=50-80' \
'ppi.hotspot_res=[A221,A222,A228,A248,A266,B140,B141,B142,B150,B151]' \
inference.ckpt_override_path=$RFD_path/models/Complex_base_ckpt.pt \
'potentials.guiding_potentials=["type:binder_ROG,weight:5,min_dist:5", "type:binder_ncontacts,weight:10", "type:interface_ncontacts,weight:5"]' \
potentials.guide_scale=0.5 \
potentials.guide_decay="linear"

# ---------------------- ProteinMPNN ----------------------

# Modified script to allow multiple target chains
../scripts/3-run_pMPNN_batch_multiple_targets.sh -i samples -o output_pMPNN -b C -t "A B" -n 1 -T 0.0001

# ---------------------- AF2 ----------------------

# Desactivate any env and activate env to run AF2-multimer (DiscobaMultimer)
conda deactivate; conda activate DiscobaMultimer

# Generate the input for DiscobaMultimer (implementation of LocalColabFold)
python ../scripts/4-generate_discoba_multimer_input_multiple_targets.py -i output_pMPNN -t input/TbTfR.fa -o discoba_multimer_designs -p input/input_structure_clean.pdb

# The latter command generated a directory. cd to it
cd discoba_multimer_designs

# Execute AF2
discoba_multimer_batch -a -i merged_MSA -c AF2_templates.config database.fasta IDs_binders.txt | tee report.log

# ---------------------- Score models ----------------------

# Deactivate environment and activate one that can permita manipulate tables and strures (here we use MultimerMapper)
conda deactivate; conda activate MultimerMapper

# cd to the project directory once more

# Analyze desings
python ../scripts/6a-analyze_binders_trf.py --af_folder discoba_multimer_designs/AF2 --rf_folder samples --output_dir analysis_plots
