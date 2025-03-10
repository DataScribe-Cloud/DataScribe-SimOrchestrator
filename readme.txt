## Please read before use:

Step 0. Place the input parameter file in the “input” folder. The input parameter file is often a "CSV" file that is created using Design Of Experiment (DOE) methods such as "LHS", "Sobol", etc. There is a small example in the "input" folder. You may use the "input.m" matlab file to create an sample-set.

Step 1. Place the simulation engine (all required files for running one instance of simulation) in the “1_setup/black_box_function” and run “Surrogate.m” in the supercomputer environment to create the simulation environment. This code will create a directory named “simulation” in the application root directory with appropriate individual folders. The current setup works for any python black-box computing engine and is tested for Fortran, Matlab, and other codes. If your code is not python, rename the appropriate lines in the Matlab code. Let me know if you need any help with this.

Note: The simulation code (black-box function) should read one set of parameters from an input file and write the output "quantities of interest" to a result text file.

Step 2. The code will also create a file named “commands.in” and will place it in the “2_running” folder. Check the example file in the folder.

Step 3. Edit "2_running/batchFile.slurm” according to your need and run it. It will execute all commands located in the “commands.in”. 

Step 4. After all of the simulations are finished, compile the individual result files generated during simulation run-time in the individual folders inside the simulation directory. 

### How to cite High-Throuput black-box function runner

# Please cite any of the below works if you like this work to help promote my work:

@article{attari2021machine,
  title={Machine Learning-Assisted High-Throughput Exploration of Interface Energy Space in Multi-Phase-Field Model with CALPHAD potential},
  author={Attari, Vahid and Arroyave, Raymundo},
  journal={Materials Theory},
  volume={6},
  pages={5},
  year={2022},
  doi={https://doi.org/10.1186/s41313-021-00038-0}
}

@article{attari2020uncertainty,
  title={Uncertainty propagation in a multiscale CALPHAD-reinforced elastochemical phase-field model},
  author={Attari, Vahid and Honarmandi, Pejman and Duong, Thien and Sauceda, Daniel J and Allaire, Douglas and Arroyave, Raymundo},
  journal={Acta Materialia},
  volume={183},
  pages={452--470},
  year={2020},
  publisher={Elsevier}
}

application structure/
│── main.py           # Entry point for the application
│── config.py         # Configuration file 
│── requirements.txt  # Dependencies list
│── README.md         # Project documentation
│── scripts/          # Batch processing scripts
│   │── batchFile.slurm  
│   │── slurm_script.sh  
│── src/  # Source code for the application
│   │── __init__.py  
│   │── black_box/  # Black-box related functions
│   │   │── __init__.py  
│   │   │── black_box_function.py  
│   │── utils/  # Helper functions
│   │   │── __init__.py  
│   │   │── create_folder_structure.py  
│   │   │── create_command_file.py  
│   │── processing/  # Data processing modules
│   │   │── __init__.py  
│   │   │── compile_data.py  
│   │   │── input_sample_generation.py  
│   │   │── samples_gen.py  
│   │── visualization/  # Data visualization
│   │   │── __init__.py  
│   │   │── data_visualization.py  
│   │   │── microstructure_visualization.py  
│── logs/  # Store logs
│   │── img_cr_log.txt  
│   │── job_log.txt  
│── input/  # Store input data
│── data/   # Processed or intermediate data
│── outputs/  # Store generated outputs
│── docs/  # Documentation

###

attari.v@C7F9W7HF-371e6c-3 1_app % python3 main.py                                     
SLURM script 'slurm_script.sh' created successfully.
#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                   # Do not propagate environment
#SBATCH --get-user-env=L                # Replicate login environment

## NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=HTP_Phasefield_LAUNCHER       # Set the job name to "JobExample3"
#SBATCH --time=0-5:00:00                         # Set the wall clock limit to 15 hours
#SBATCH --ntasks=2000                            # Request 500 tasks
#SBATCH --ntasks-per-node=40                     # Request 20 tasks/cores per node
#SBATCH --cpus-per-task=1                        # Request 1 CPU per task
###SBATCH --mem-per-cpu=3000                     # Commented out memory per CPU
#SBATCH --mem=40000M                             # Request 40GB per node
#SBATCH --output=code_output%J.txt               # Send stdout/err to "Example3Out.[jobID]"
#SBATCH --account=XXXXXXXXXXXX                   # Specify account number

ml GCCcore/12.3.0

tamulauncher --remove-logs --commands-pernode=40 commands.in
tamulauncher commands.in

Welcome to Task Manager!

What would you like to do?
1. Sample generation
2. Read existing sample file
3. Quit
Enter your choice (1:SampleGeneration/2:Read/3:Quit): 2
Existing sample read...
           eps       tav         K        Te     alpha        nu  Heat_diffusivity
0     0.023334  0.000196  1.599945  1.077113  0.899420  0.035152          0.854491
1     0.011431  0.000193  1.008234  0.844074  1.127658  0.009141          1.087716
2     0.021329  0.000176  1.894478  0.969898  1.095646  0.016931          0.852770
3     0.020397  0.000170  1.392461  0.835690  1.120068  0.001203          1.029000
4     0.029577  0.000165  1.325258  1.164922  0.942478  0.005756          1.080566
...        ...       ...       ...       ...       ...       ...               ...
4995  0.009217  0.000122  0.822715  0.949437  0.988719  0.001298          1.044566
4996  0.007294  0.000142  1.516017  1.008614  1.163013  0.036333          1.119603
4997  0.009584  0.000157  0.531135  1.001451  1.118030  0.038402          1.060784
4998  0.016651  0.000129  1.061001  1.050996  0.805083  0.030338          1.013816
4999  0.019264  0.000102  0.848054  0.913693  1.194965  0.001528          0.934004

[5000 rows x 7 columns]

 **************** TAMULAUNCHER - Compiling app **************** 
 ************************************************************** 
  

Make command output:
rm -f kobayashi *.o *.mod *.txt *.out *.mp4
gfortran -ffree-line-length-none -O2 -Wall -o kobayashi Kobayashi.f90

Make command executed successfully!


 **************** Creating the folder structure ***************
 ************************************************************** 
  
File samples_lhs.csv is retrieved and samples will be used for calculations.

 **************** detect the application files ************ 
 **************** detect the application files ************ 
['Kobayashi.f90', 'mod_geometry.mod', '.DS_Store', 'makefile', 'wrt_opts.mod', 'heat_eqn.mod', 'pfm.mod', 'mod_electromigration.mod', 'math_opts.mod', 'kobayashi', 'mts_pars.mod']

Executable file(s) were detected as: kobayashi

 **************** TAMULAUNCHER ENV. BUILDER ************** 
 ********************************************************* 
 ********************************************************* 
 **** Starting parameter-set                         :      1
 **** Last parameter-set                             :      5000
total_iterations: 5000
Processing:   0%|                                                                                                                                                                                                        | 0/5000 [00:00<?, ?iteration/s]
Would you like to delete old simulations directory?
1. Delete Directory
2. Don't Delete Directory
Enter your choice (1:Yes/2:No): 2
Exiting Task Manager.
Processing: 100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 5000/5000 [01:31<00:00, 54.71iteration/s]
Loop completed

 **************** TAMULAUNCHER Command FL BUILDER ************* 
 ************************************************************** 
  
Commands written to commands.in
Welcome to Task Manager!

What would you like to do?
0. Perform Local run
1. Perform Terra run
2. Perform Grace run
3. Perform Faster run
4. Quit
Enter your choice (0:Local/1:Terra/2:Grace/3:Faster/4:Quit): ???