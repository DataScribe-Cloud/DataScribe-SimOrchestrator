#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                   # Do not propagate environment
#SBATCH --get-user-env=L                # Replicate login environment

## NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=HTP_Phasefield_LAUNCHER       # Set the job name to "JobExample3"
#SBATCH --time=0-5:00:00                        # Set the wall clock limit to 15 hours
#SBATCH --ntasks=2000                            # Request 500 tasks
#SBATCH --ntasks-per-node=40                     # Request 20 tasks/cores per node
#SBATCH --cpus-per-task=1                        # Request 1 CPU per task
###SBATCH --mem-per-cpu=3000                     # Commented out memory per CPU
#SBATCH --mem=40000M                             # Request 40GB per node
#SBATCH --output=code_output%J.txt               # Send stdout/err to "Example3Out.[jobID]"
#SBATCH --account=132794105877                   # Specify account number

ml GCCcore/12.3.0

tamulauncher --remove-logs --commands-pernode=40 commands.in
tamulauncher commands.in
