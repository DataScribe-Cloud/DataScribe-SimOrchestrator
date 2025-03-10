# code
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import os
import sys
import re
from joblib import Parallel, delayed

import logging

# Configure logging
logging.basicConfig(filename='job_log.txt', level=logging.INFO, format='%(asctime)s - %(message)s')

APP_FOLDER = os.getcwd()

def read_vahid_file(filename):
    file_obj = open(filename,'r')
    zone_dict = {}
    array_list = []
    for line in file_obj:
        if 'variables' in line:
            re_finds_obj = re.findall(r'["](.*?)["]',line)
            headers = [find for find in re_finds_obj]
            
        elif 'ZONE' in line:
            # pass
            re_finds_obj = re.findall(r'[I][=][\s]*[0-9]+',line)
            I_val = int(re_finds_obj[0].split('=')[1].strip())
            #print('re_finds_obj',re_finds_obj)
            re_finds_obj = re.findall(r'[J][=][\s]*[0-9]+',line)
            J_val = int(re_finds_obj[0].split('=')[1].strip())
            zone_dict['I'] = I_val
            zone_dict['J'] = J_val
        elif line:
            line = line.replace('D','E')
            fix_E_points = re.findall('[\d][-][\d]',line)
            for E_point in fix_E_points:
                new_E_point = E_point.replace('-','E-')
#                print(new_E_point)
                line = line.replace(E_point,new_E_point)
#            print(line)
            #line = re.sub('^[\s]+','',line)
            #line = re.sub('[\s]+',' ', line)
            # print(line)
            
            line = np.fromstring(line.strip('[]'), dtype=float, sep = ' ')
            
            array_list.append(line)
            #array_list.append(np.fromstring(line,sep = ' '))
        else:
            pass
    #variables =  "IG"  "JG"  "PHIMAX"  "C1"  "PHI1_2"
    #file_obj.readlines()
    
    
                    
    file_obj.close()
    df_out = pd.DataFrame(array_list,columns=headers)
    return df_out, zone_dict

def read_samples():
	
	samples = pd.read_csv('../input/samples_lhs.csv')
	
	return samples


def read_comp(fl_name,id):
			
	#print('folder ',id)
	
	try:
		fl = os.path.join(APP_FOLDER,'..','simulations',str(id),'results',fl_name)		
		df, zone = read_vahid_file(fl)
		#print(df)
		#print(zone)
		#print(zone['I'])
		#print(zone['J'])
		
		comp = df['PHI'].tolist()
		#print(comp)
		N = int(zone['I']) #int(np.sqrt(np.size(comp)))
		M = int(zone['J']) #int(np.sqrt(np.size(comp)))
		img = np.flip (np.transpose( np.reshape(np.array(comp), (M, N) ) ) )
	except Exception as error:
		N = 401
		M = 101
		img = np.flip( np.transpose( np.zeros([N, M]) ) )
		print("An exception occurred:", error) 
		
	return img

def plot_images(id):
	
	logging.info(f"Processing ID: {id}")
	
	fl_name = 'PhiTemp000005000.plt'
		
	#print('folder ',id)
	
	try:
		img = read_comp(fl_name,id)
		plt.imshow(img,vmin=0.17, vmax=0.86 )
		# set axis spines (the box around the plot) to be invisible
		ax = plt.gca()
		ax.axes.get_xaxis().set_ticks([])
		ax.axes.get_yaxis().set_ticks([])
		plt.text(6, 18, str(id), bbox=dict(fill='white', edgecolor='blue', linewidth=1))
		plt.savefig('../data/mics/mic_'+fl_name[:-4]+'_'+str(id)+'.jpg', dpi=300, bbox_inches = 'tight', pad_inches = 0) #, pad_inches = 0
	except Exception as error:
		print("An exception occurred:", error) 


if __name__ == "__main__":

	APP_FOLDER = os.getcwd()
	for dirs in os.walk(os.path.join(APP_FOLDER,'..','simulations')):
		totalDir=0
		for directories in dirs:
			totalDir += 1
	        
	print(totalDir)
	
	## read input samples
	samples = read_samples()
	## read input samples

	## ADD NUMBER OF FOLDERS HERE
	Folders = np.size(samples,0) #int(input("Enter Number of Simulations (Positive Integer ()): "))
	
	print(Folders)
	
	dest_dir = '../data/mics'

	try:
		if os.path.exists(dest_dir):
			shutil.rmtree(dest_dir)
			print('replacing the src folder - ',str(id) )
	except Exception as error:
		print("An exception occurred:", error) 
            
	try:
		#os.system('mkdir ../data/mics')
		os.makedirs('../data/mics', exist_ok=True)
	except:
		print('mics folder already exists...')

		
	Parallel(n_jobs=6)(delayed(plot_images)(id) for id in range(0,5000))

		
	
	
