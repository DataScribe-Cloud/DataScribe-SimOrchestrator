import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
from joblib import Parallel, delayed
import logging
import shutil
import time

# Configure logging
logging.basicConfig(filename='job_log.txt', level=logging.INFO, format='%(asctime)s - %(message)s')

APP_FOLDER = os.getcwd()

def read_vahid_file(filename):
    with open(filename, 'r') as file_obj:
        zone_dict = {}
        array_list = []
        headers = []
        
        for line in file_obj:
            if 'variables' in line:
                headers = re.findall(r'["](.*?)["]', line)
            elif 'ZONE' in line:
                I_val = int(re.search(r'I\s*=\s*(\d+)', line).group(1))
                J_val = int(re.search(r'J\s*=\s*(\d+)', line).group(1))
                zone_dict['I'] = I_val
                zone_dict['J'] = J_val
            elif line.strip():
                line = line.replace('D', 'E')
                line = re.sub(r'(\d)-(\d)', r'\1E-\2', line)
                line_data = np.fromstring(line.strip('[]'), dtype=float, sep=' ')
                array_list.append(line_data)
        
    df_out = pd.DataFrame(array_list, columns=headers)
    return df_out, zone_dict

def read_samples():
    return pd.read_csv('../input/samples_lhs.csv')

def read_comp(fl_name, id, retry_count=5, delay=2):
    logging.info(f"Reading computation for folder {id}")
    fl = os.path.join(APP_FOLDER, '..', 'simulations', str(id), 'results', fl_name)
    
    for attempt in range(retry_count):
        try:
            logging.info(f"Checking file path: {fl}")
            if not os.path.isfile(fl):
                raise FileNotFoundError(f"File not found: {fl}")
            
            df, zone = read_vahid_file(fl)
            comp = df['PHI'].to_numpy()
            
            N, M = zone['I'], zone['J']
            #if comp.size != N * M:
            #    print(comp)
            #    raise ValueError(f"Cannot reshape array of size {comp.size} into shape ({M},{N})")
            
            img = np.flip(np.transpose(comp.reshape((M, N))))
            return img
        
        except (FileNotFoundError, ValueError) as error:
            #logging.error(f"Attempt {attempt + 1} - Error in folder {id}: {error}")
            time.sleep(delay)
    
    logging.error(f"Failed to read computation for folder {id} after {retry_count} attempts")
    N, M = 101, 401
    return np.flip(np.transpose(np.zeros((N, M))))

def plot_images(id):
    logging.info(f"Processing ID: {id}")
    fl_name = 'PhiTemp000005000.plt'
    try:
        img = read_comp(fl_name, id)
        plt.imshow(img, vmin=0.17, vmax=0.86)
        ax = plt.gca()
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])
        plt.text(6, 18, str(id), bbox=dict(fill='white', edgecolor='blue', linewidth=1))
        plt.savefig(f'../data/mics/mic_{fl_name[:-4]}_{id}.jpg', dpi=300, bbox_inches='tight', pad_inches=0)
        plt.close()
    except Exception as error:
        logging.error(f"An error occurred while plotting image for ID {id}: {error}")

if __name__ == "__main__":
    APP_FOLDER = os.getcwd()
    
    #for dirs in os.walk(os.path.join(APP_FOLDER, '..', 'simulations')):
    #    totalDir = sum([len(directories) for directories in dirs])
    
    #print(totalDir)
    
    samples = read_samples()
    Folders = len(samples)
    print(Folders)
    
    dest_dir = '../data/mics'
    if os.path.exists(dest_dir):
        shutil.rmtree(dest_dir)
        print(f'Replacing the src folder - {dest_dir}')
    
    os.makedirs(dest_dir, exist_ok=True)
    
    start = 0
    end   = 5000
    print('start:',start)
    print('end  :',end)
    
    Parallel(n_jobs=8)(delayed(plot_images)(id) for id in range(start, end))
    
