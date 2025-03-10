#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 17:39:48 2023

@author: attari.v
"""

## compile_app

import os
import subprocess

def run_make_command(srcdir):


    print(' ');
    print(' **************** TAMULAUNCHER - Compiling app **************** ');
    print(' ************************************************************** \n');    

    try:

        new_path = os.path.join(srcdir,'black_box_function','src')
        os.chdir(new_path)
        
        # Run the make command
        completed_process = subprocess.run(["make clean; make;"], capture_output=True, text=True, check=True, shell=True)

        # Print the captured stdout and indicate success
        print("")
        print("Make command output:")
        print(completed_process.stdout)
        print("Make command executed successfully!")
        print("")

        os.chdir(srcdir)

    except subprocess.CalledProcessError as e:
        # Print the captured stderr and indicate failure
        print("Make command error:")
        print(e.stderr)
        print("Make command failed!")
