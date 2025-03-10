#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 22:16:04 2023

@author: attari.v ad
"""

import os
import subprocess
#from input_sample_generation import _lhs_init
from .input_sample_generation import _lhs_init
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def samples_gen():
    print("Welcome to Task Manager!")
    
    while True:
        print("\nWhat would you like to do?")
        print("1. Sample generation")
        print("2. Read existing sample file")
        print("3. Quit")
        
        choice = input("Enter your choice (1:SampleGeneration/2:Read/3:Quit): ")
        
        if choice == '1':
            parameters = perform_task_a()
            break  # exit the loop
        elif choice == '2':
            parameters = perform_task_b()
            break  # exit the loop
        elif choice == '3':
            print("Exiting Task Manager.")
            break  # exit the loop
        else:
            print("Invalid choice. Please select a valid option.")

    return parameters            


def perform_task_a():
    print("Sample generation...")

    try:

        num_samples = input("Enter Number of Desired Samples (Positive Integer (>1)): ") 

        # ## LHS design     
        par_names = ['eps','tav','K','Te','alpha','nu','Heat_diffusivity']
        bounds    = np.log10([(7e-3,3e-2),(1e-4,2e-4),(0.5,2.2),(0.8,1.2),(0.8,1.2),(0.001,0.05),(0.85,1.15)])
        num_samples   = int(num_samples)
    
        parameters = _lhs_init(par_names, bounds, num_samples)
        parameters = 10**(parameters); # or exp
    
        print(parameters)
    
        fig, ax = plt.subplots(figsize = (9, 6))
        # ax.scatter(parameters.Ls,parameters.Le, alpha=0.7, edgecolors="k")
        # Set logarithmic scale on the both variables
        ax.set_xscale("log")
        ax.set_yscale("log");
        
        os.makedirs('../input', exist_ok=True)
    
        parameters.to_csv('../input/samples_lhs.csv',index=False)   

    except Exception as error:
        print("An exception occurred:", error) # An exception occurred: division by zero
		
    return parameters


def perform_task_b():
    print("Existing sample read...")
    
    parameters = pd.read_csv('../input/samples_lhs.csv')
    
    print(parameters)
    
    return parameters


def perform_task_c():
    print("Task HTP-TAMULAUNCHER run on Faster.")


def run_command_in_background(command, log_file):
    with open(log_file, 'w') as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    
    return process
