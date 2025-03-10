import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import glob
import scipy.stats

import csv

plt.rcParams.update({'font.size': 15})

# %%

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
			#print(re_finds_obj)
			#re_finds_obj = re.findall(r'[J][=][\s]*[0-9]+',line)
			J_val = int(re_finds_obj[0].split('=')[1].strip())
			zone_dict['I'] = I_val
			zone_dict['J'] = J_val
		elif line:
			line = line.replace('D','E')
			#line = re.sub('^[\s]+','',line)
			#line = re.sub('[\s]+',' ', line)
			line = np.fromstring(line.strip('[]'), sep = ' ')
			array_list.append(line)
			#array_list.append(np.fromstring(line,sep = ' '))
		else:
			pass
	#variables =  "IG"  "JG"  "PHIMAX"  "C1"  "PHI1_2"
	#file_obj.readlines()
	
	file_obj.close()
	df_out = pd.DataFrame(array_list,columns=headers)
	
	return df_out, zone_dict    
	
# %%

## Read the result files for three time iterations, calculate the max. and min compositions
## and save into csv files... Find non-spinodal cases

cols = ['composition']

fl_names = ['phi000010000.plt']


os.makedirs('images10000',exist_ok=True)

for fl in fl_names:

	print()
	print(fl)
	
	diff = []
	max_comp = []
	min_comp = []
	for i in range(0,8):
	
		if np.remainder(i, 1000)==0:
			print(i)
			
		# file path to read 
		fl_path = os.path.join('simulations',str(i),'input','var.csv')
		with open(fl_path, newline='') as csvfile:
			spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
			for row in spamreader:
				print(', '.join(row))
				

		# file path to read 
		fl_path = os.path.join('simulations',str(i),'microstructure',fl)
		#print(fl_path)
		# read the file and identify the non-spinodal cases
		df, zone_dict = read_vahid_file(fl_path)
		
		array = np.array(df.phi1.values)
		N = int(np.sqrt(np.size(array)))
		img = np.reshape(np.array(array), (N, N) )
		plt.imshow(img)
		plt.axis('off')
		plt.savefig('images10000/img_'+str(fl[:-4])+'_'+str(i)+'.png', bbox_inches = 'tight',pad_inches = 0)
		plt.show()
		
		#max_c = max(array)
		#min_c = min(array)
		#difference = (max_c - min_c)
		