#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                 			  #Do not propagate environment
#SBATCH --get-user-env=L                          #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=HTP_Phasefield_LAUNCHER       #Set the job name to "JobExample3"
#SBATCH --time=2-00:00:00                          #Set the wall clock limit to 1 Day and 12hr
#SBATCH --ntasks=500                 			  #Request 8 tasks
#SBATCH --ntasks-per-node=40        			  #Request 2 tasks/cores per node
#SBATCH --cpus-per-task=1				
###SBATCH --mem-per-cpu=3000	
#SBATCH --mem=40000M            				  #Request 4096MB (4GB) per node 
#SBATCH --output=code_output%J.txt  			  #Send stdout/err to "Example3Out.[jobID]"
#SBACTH --account=132794105877

#load_pandas
#load_numpy
#load_mls
#ml intel/2021a

ml GCCcore/12.3.0

tamulauncher --remove-logs commands.in
tamulauncher commands.in


