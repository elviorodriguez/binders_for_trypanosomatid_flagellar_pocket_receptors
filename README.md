# Binders for trypanosomatid flagellar pocket receptors
This repository contains the code to reproduce the designs of binders targeted against flagellar pocket receptors from trypanosomatids. It also contains the analysis of the resulting designs (3D scatter plots and 2D scatter plots).

## Design directories
The directories ```hphbr_1```, ```hphbr_2```, ```isg65``` and ```tfr``` contain the code and files needed to reproduce the designs from the manuscript. Inside each directory you can find:

  - ```design_commands.sh```: Contains the code that need to be executed in the command line (we used bash inside Ubuntu 24.04). The code is split in sections that correspond to the different steps of the design process (RFdiffusion, ProteinMPNN, AF2-multimer and Scoring).
  - ```input```: Directory containing the input files for each of the steps (PDBs and FASTAs)
  - ```analysis_plots```: Contains the analysis of the resulting designs. ```binder_analysis_report.html``` is an HTML report that summarizes all the plots interactively.

## Scripts directory
The directory ```scripts``` contains all the scripts that are called by ```design_commands.sh```. Some of them uses specific conda environments to run each of the different softwares (RFdiffusion, ProteinMPNN and AF2-multimer). If you want to reproduce the results or use the scripts with your own data, you will need to first create and then activate this environments in the order described in ```design_commands.sh```. Have a look at these repos to see how to install them:

  - RFdiffusion: https://github.com/RosettaCommons/RFdiffusion
  - ProteinMPNN: https://github.com/dauparas/ProteinMPNN
  - AF2-multimer (using DiscobaMultimer): https://github.com/elviorodriguez/DiscobaMultimer

## R analysis code
The directory ```metrics_analysis_R``` contains the code to reproduce the 2D scatter plot from the manuscript.
