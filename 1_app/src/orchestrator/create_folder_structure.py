import os
import pandas as pd
import numpy as np
import shutil
import subprocess
from tqdm import tqdm
import time

def create_folder_structure(samples):
    print(' ')
    print(' **************** Creating the folder structure ***************')
    print(' ************************************************************** \n')    

    srcdir = os.getcwd()
    wrkdir = os.path.join(srcdir, '../..')  # Change the number

    inputdir = os.path.join(wrkdir, 'input')
    datadir = os.path.join(wrkdir, 'data')

    try:
        for file in os.listdir(inputdir):
            if file.endswith(".csv"):
                csv_file = file
        filename = os.path.join(inputdir, csv_file)
        samples = pd.read_csv(filename)
        
        print(f'File {csv_file} is retrieved and samples will be used for calculations.')
        print()
        
    except Exception as error:
        print(csv_file)  # file list 
        print("An exception occurred:", error)  # An exception occurred: division by zero
        print("An exception occurred:", type(error).__name__)  # An exception occurred: ZeroDivisionError
    
    print(' **************** detect the application files ************ ')
    print(' **************** detect the application files ************ ')

    try:
        src_path = os.path.join(srcdir, 'black_box_function', 'src')
        app_fls = os.listdir(src_path)
        print(app_fls)
        
        executables = []
        for fl in app_fls:
            fl_path = os.path.join(src_path, fl)
            if os.path.isfile(fl_path) and os.access(fl_path, os.X_OK):
                executables.append(fl)
                print()
                print('Executable file(s) were detected as:', fl)
                print()
                    
    except Exception as error:
        print("An exception occurred:", error)  # An exception occurred: division by zero
        print("An exception occurred:", type(error).__name__)  # An exception occurred: ZeroDivisionError
    
    startRow = 1  # Folder 1
    endRow = np.size(samples, 0)  # iFolder "end"

    print(' **************** TAMULAUNCHER ENV. BUILDER ************** ')
    print(' ********************************************************* ')
    print(' ********************************************************* ')
    print(' **** Starting parameter-set                         :     ', startRow)
    print(' **** Last parameter-set                             :     ', endRow)

    try:
        os.makedirs(os.path.join(wrkdir, 'simulations'), exist_ok=True)
    except Exception as error:
        print("An exception occurred:", error)  # An exception occurred: division by zero

    total_iterations = np.size(samples, 0)
    print('total_iterations:',total_iterations)

    progress_bar = tqdm(total=total_iterations, desc="Processing", unit="iteration")

    simulationsdir = os.path.join(wrkdir, 'simulations')
    try:
        if os.path.exists(simulationsdir):
            dir_deletion_request()
    except Exception as error:
        print("An exception occurred:", error)  # An exception occurred: division by zero
            
    os.makedirs(simulationsdir, exist_ok=True)
    
    for i in range(0, total_iterations):
        try:    
            os.chdir(simulationsdir)
            os.makedirs(str(i), exist_ok=True)
            os.makedirs(os.path.join(str(i), 'input'), exist_ok=True) 
        except Exception as error:
            print("An exception occurred:", error)  # An exception occurred: division by zero
            
        inputdir = os.path.join(simulationsdir, str(i), 'input')
        
        sample = samples.iloc[i, :]
        sample.transpose().to_csv(os.path.join(inputdir, 'var.csv'), index=False, header=False, sep=',')
        #sample.T.to_csv(os.path.join(inputdir, 'var.csv'), index=False, sep=',')
        #sample.T.reset_index(drop=True).to_csv(os.path.join(inputdir, 'var.csv'), index=False, sep=',')

        try:
            src_dir = os.path.join(srcdir, 'black_box_function', 'src')
            dest_dir = os.path.join(simulationsdir, str(i), 'src')
            if os.path.exists(dest_dir):
                shutil.rmtree(dest_dir)
                print('replacing the src folder - ', str(i))
            shutil.copytree(src_dir, dest_dir)
        except Exception as error:
            print("An exception occurred:", error)  # An exception occurred: division by zero

        progress_bar.update(1)        
    
    try:
        src_fl = os.path.join(srcdir, 'batchFile.slurm')
        dest_fl = os.path.join(simulationsdir, 'batchFile.slurm')
        shutil.copy(src_fl, dest_fl)
    except Exception as error:
        print("An exception occurred:", error)  # An exception occurred: division by zero

    progress_bar.close() 
    print("Loop completed")

def dir_deletion_request():
    while True:
        print("\nWould you like to delete old simulations directory?")
        print("1. Delete Directory")
        print("2. Don't Delete Directory")
        
        choice = input("Enter your choice (1:Yes/2:No): ")
        
        if choice == '1':
            perform_task_a()
            break
        elif choice == '2':
            print("Exiting Task Manager.")
            break
        else:
            print("Invalid choice. Please select a valid option.")

def perform_task_a():
    print("Perform directory deletion...")
    
    print(os.getcwd())

    try:
        simulationsdir = '../../simulations'

        if os.path.exists(simulationsdir):
            shutil.rmtree(simulationsdir)
            print('OLD simulations directory deleted - ')

        os.chdir('../../1_app')
        

    except subprocess.CalledProcessError as e:
        print(e.stderr)

def perform_task_b():
    print("Task HTP-TAMULAUNCHER run on Terra.")
    