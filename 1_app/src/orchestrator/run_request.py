#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 22:16:04 2023

@author: attari.v
"""

import os
import subprocess

def run_request(num_samples):
    print("Welcome to Task Manager!")
    
    while True:
        print("\nWhat would you like to do?")
        print("0. Perform Local run")
        print("1. Perform Terra run")
        print("2. Perform Grace run")
        print("3. Perform Faster run")
        print("4. Quit")
        
        choice = input("Enter your choice (0:Local/1:Terra/2:Grace/3:Faster/4:Quit): ")
        
        if choice == '0':
            perform_task_local(num_samples)
            break  # exit the loop
        elif choice == '1':
            perform_task_a()
            break  # exit the loop
        elif choice == '2':
            perform_task_b()
            break  # exit the loop
        elif choice == '3':
            perform_task_c()
            break  # exit the loop
        elif choice == '4':
            print("Exiting Task Manager.")
            break  # exit the loop
        else:
            print("Invalid choice. Please select a valid option.")

def perform_task_local(num_samples):
    print("Task HTP-serial run on Local Computer.")
    
    try:
        
        #os.chdir('../simulations')
        print(os.getcwd())
        
        for idx in range(num_samples):
            
            print('running simulation '+str(idx)+'/src')
                        
            os.chdir(str(idx)+'/src')
            log_file = "batchFile-output.log"
            command  = "./kobayashi;"
            background_process = run_command_in_background(command, log_file)
            os.chdir('../..')
        
        #with open(log_file, 'w') as log:
			# Run the make command
            #completed_process = subprocess.Popen(["sbatch batchFile.slurm"], stdout=log, stderr=subprocess.STDOUT, text=True, shell=True)
        
        #completed_process = subprocess.run(["tamulauncher batchFile.slurm;"], capture_output=True, text=True, check=True, shell=True)
        #completed_process = subprocess.Popen(["tamulauncher batchFile.slurm; &"], shell=True)
        
        
        # Print the captured stdout and indicate success
        print(" ")
        #print("Output:")
        #print(completed_process.stdout)
        print("Local calculations Executed Successfully in the Background!")
        print(" ")
        
        os.chdir('../1_app')
        
    except subprocess.CalledProcessError as e:
        # Print the captured stderr and indicate failure
        print("Make command error:")
        print(e.stderr)
        print("Make command failed!")
        
def perform_task_a():
    print("Task HTP-serial run on Terra.")

    try:

        os.chdir('../simulations')
		
        if os.path.exists(".tamulauncher-log"):
            os.rmdir([".tamulauncher-log"])
		
        log_file = "batchFile-output.log"
        command  = "sbatch batchFile.slurm;"
        background_process = run_command_in_background(command, log_file)

        #with open(log_file, 'w') as log:
			# Run the make command
            #completed_process = subprocess.Popen(["sbatch batchFile.slurm"], stdout=log, stderr=subprocess.STDOUT, text=True, shell=True)

        #completed_process = subprocess.run(["tamulauncher batchFile.slurm;"], capture_output=True, text=True, check=True, shell=True)
        #completed_process = subprocess.Popen(["tamulauncher batchFile.slurm; &"], shell=True)

		
        # Print the captured stdout and indicate success
        print(" ")
        #print("Output:")
        #print(completed_process.stdout)
        print("Tamulauncher Command Executed Successfully in the Background!")
        print(" ")

        os.chdir('../1_app')
        
    except subprocess.CalledProcessError as e:
        # Print the captured stderr and indicate failure
        print("Make command error:")
        print(e.stderr)
        print("Make command failed!")



def perform_task_b():
    print("Task HTP-TAMULAUNCHER run on Grace.")
    
    print("Creating batchSlurm.sbatch file for Grace. Edit the file parameters on 'run_request.py'. ")
    create_slurm_file(filename="batchFile.slurm")
    
    try:
        
        os.chdir('../simulations')
        
        if os.path.exists(".tamulauncher-log"):
            os.rmdir([".tamulauncher-log"])
            
        log_file = "batchFile-output.log"
        command  = "sbatch batchFile.slurm;"
        background_process = run_command_in_background(command, log_file)
                
        # Print the captured stdout and indicate success
        print(" ")
        #print("Output:")
        #print(completed_process.stdout)
        print("Tamulauncher Command Executed Successfully in the Background!")
        print(" ")
        
        os.chdir('../1_app')
        
    except subprocess.CalledProcessError as e:
        # Print the captured stderr and indicate failure
        print("Make command error:")
        print(e.stderr)
        print("Make command failed!")



def perform_task_c():
    print("Task HTP-TAMULAUNCHER run on Faster.")


def run_command_in_background(command, log_file):
    with open(log_file, 'w') as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, shell=True, text=True)
    
    return process



def create_slurm_file(filename="slurm_script.sh"):
    slurm_content = """#!/bin/bash
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
"""

    with open(filename, "w") as file:
        file.write(slurm_content)
        
    print(f"SLURM script '{filename}' created successfully.")


    with open(filename, "r") as file:
        content = file.read()
        print(content)
    
# Call the function to create the SLURM script
create_slurm_file()