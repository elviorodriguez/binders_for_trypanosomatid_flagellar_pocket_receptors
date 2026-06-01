# cd to the directory to the project and work from there

# Set your own paths:
RFD_path=/home/esteban/RFdiffusion

# ---------------------- RFdiffusion ----------------------

# Activate RFdiff environment
conda activate SE3nv

# Run RFdiff
$RFD_path/scripts/run_inference.py 'inference.input_pdb=/home/elvio/Desktop/hphbr/input/hphbr_clean.pdb' 'contigmap.contigs=[D38-297/0 10-35/C300-305/5-30/C275-280/10-35]' 'inference.num_designs=1000' 'contigmap.length=70-90' 'contigmap.inpaint_seq=[C276,C277,C279,C301,C302,C303,C304]' 'ppi.hotspot_res=[D70,D71,D72,D73,D74,D75,D78,D79,D81,D82,D85,D86,D251,D252,D255,D256,D258,D259,D262]' inference.ckpt_override_path=$RFD_path/models/Complex_base_ckpt.pt 'potentials.guiding_potentials=["type:binder_ROG,weight:5,min_dist:5","type:binder_ncontacts,weight:5", "type:interface_ncontacts,weight:10"]' potentials.guide_scale=0.5 potentials.guide_decay="linear"

# ---------------------- ProteinMPNN ----------------------

# Run ProteinMPNN (ProteinMPNN environment is activated within run_pMPNN_batch.sh)
../scripts/3-run_pMPNN_batch.sh -i samples -o output_pMPNN -b "C" -t "D"

# ---------------------- AF2 ----------------------

# Desactivate any env and activate env to run AF2-multimer (DiscobaMultimer)
conda deactivate; conda activate DiscobaMultimer

# Generate the input for DiscobaMultimer (implementation of LocalColabFold)
python ../scripts/4-generate_discoba_multimer_input.py -i output_pMPNN -t input/hphbr.fa -o discoba_multimer_designs -p input/hphbr_clean.pdb

# The latter command generated a directory. cd to it
cd discoba_multimer_designs

# Execute AF2
discoba_multimer_batch -a -i merged_MSA -c AF2_templates.config database.fasta IDs_binders.txt | tee report.log


# ---------------------- Score models ----------------------

# Deactivate environment and activate one that can permita manipulate tables and PDB structures with BioPython (here we use MultimerMapper)
conda deactivate; conda activate MultimerMapper

# Analyze desings
python ../scripts/6b-analyze_binders.py --af_folder discoba_multimer_designs/AF2 --rf_folder samples --output_dir analysis_plots --binder_chain_rf "C" --binder_chain_af "A"
