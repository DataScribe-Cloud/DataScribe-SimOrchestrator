#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 12:27:57 2023

@author: attari.v
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Ensure src is in the Python path
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))

# Import functions from the newly created orchestrator module
from orchestrator.samples_gen import samples_gen
from orchestrator.create_folder_structure import create_folder_structure
from orchestrator.create_command_file import create_command_file
from orchestrator.run_request import run_request
from orchestrator.compile_app import run_make_command

# Define directory paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))  # Root directory of the project
SRC_DIR = os.path.join(BASE_DIR, "src")
INPUT_DIR = os.path.join(BASE_DIR, "input")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
LOG_DIR = os.path.join(BASE_DIR, "logs")

def main():
    """Main function to run the simulation orchestration"""

    print("Initializing Simulation Orchestrator...")

    # Step 1: Generate parameters
    parameters = samples_gen()
    
    # Step 2: Compile application
    print("Compiling application...")
    run_make_command(SRC_DIR)

    # Step 3: Create necessary folder structure
    create_folder_structure(parameters)

    # Step 4: Create command file for batch processing
    create_command_file(fl_name=os.path.join(INPUT_DIR, 'commands.ini'),
                        samples=parameters, 
                        app_flname='kobayashi')

    # Step 5: Submit run request
    run_request(parameters.shape[0])

    print("Simulation Orchestrator has finished execution.")
    
if __name__ == "__main__":
    main()