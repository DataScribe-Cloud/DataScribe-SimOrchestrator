#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 12:54:49 2023

@author: attari.v
"""


import os



def create_command_file(fl_name,samples,app_flname):
    
    print(' ');
    print(' **************** TAMULAUNCHER Command FL BUILDER ************* ');
    print(' ************************************************************** \n');    

    try:

        fl_name = 'commands.in'  # File name
        with open(fl_name, 'w') as f:  # Open the file for writing
            for i in range(0, len(samples)):
                # Write to the file using the file object 'fl'
                f.write('cd ../simulations/'+str(i)+'/src; ./'+str(app_flname)+'; \n')
        
        print("Commands written to", fl_name)
    except Exception as error:
        print("An exception occurred:", error) # An exception occurred: division by zero
        print("An exception occurred:", type(error).__name__) # An exception occurred: ZeroDivisionError
