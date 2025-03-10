import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import os
import sys
import re
from tqdm import tqdm
import time

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
            fix_E_points = re.findall('[\d][-][\d]',line)
            for E_point in fix_E_points:
                new_E_point = E_point.replace('-','E-')
#                print(new_E_point)
                line = line.replace(E_point,new_E_point)
#            print(line)
            #line = re.sub('^[\s]+','',line)
            #line = re.sub('[\s]+',' ', line)
            # print(line)
            
            line = np.fromstring(line.strip('[]'), dtype=float, sep = ',')
            
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
	
	samples = pd.read_csv('../data/samples.csv')
	
	return samples


def plot_data(data):
	
	log_xcolumns = ["Ls","Le","W","eps"]
	log_ycolumns = ["features1"]

	lables_xcolumns = ["Interface Mobility","Reaction Rate","Barrier Height of Transformation","Gradient Energy Coef."]
	lables_ycolumns = ["features1: Thickness (m/s)"]

	j=0
	fig, ax = plt.subplots(1,4,sharey=True, figsize = (20,4))
	for idx_x in log_xcolumns:
		for idx_y in log_ycolumns:
			print(idx_x)
			ax[j].scatter(data[idx_x],data[idx_y], alpha=0.7, edgecolors="k", c=data.index)
			ax[j].set_xscale("log");
			#ax[j].set_yscale("log");
			ax[j].set_xlabel(lables_xcolumns[j])
			ax[j].set_ylabel(idx_y+' :Electrode Thickness (m)')
			#ax.fontsize = 15
			
			j = j+1

	fig.savefig('../data/data_'+str(log_ycolumns)+'.jpg',dpi=200)

	log_xcolumns = ["Ls","Le","W","eps"]
	log_ycolumns = ["features2"]

	j=0
	fig, ax = plt.subplots(1,4,sharey=True, figsize = (20,4))
	for idx_x in log_xcolumns:
		for idx_y in log_ycolumns:
			print(idx_x)
			ax[j].scatter(data[idx_x],data[idx_y], alpha=0.7, edgecolors="k", c=data.index)
			ax[j].set_xscale("log");
			#ax[j].set_yscale("log");
			ax[j].set_xlabel(lables_xcolumns[j])
			ax[j].set_ylabel(idx_y+str(' :Interface Velocity (m/s)'))
			#ax.fontsize = 15
			
			j = j+1

	plt.tight_layout()		
	fig.savefig('../data/data_'+str(log_ycolumns)+'.jpg',dpi=200)



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

	print(APP_FOLDER)
	symbol= '//'
	variables =  ["itimes" , "TIME" , "phitot" , "ZSOLTHK" , "thickness" , "tip_vel" , "tip_vel_analytic" , "i_n_pfm" , "i_n_pfm/i0"]

    # Number of iterations in your loop
	total_iterations = Folders
	print(total_iterations)
	# Create a tqdm instance
	progress_bar = tqdm(total=total_iterations, desc="Processing", unit="iteration")

	combined_last_rows = pd.DataFrame(columns = variables)
	log       = []
	features1 = []
	features2 = []
	features3 = []
	for x in range(int(Folders)):
		try:
			fl = os.path.join(APP_FOLDER,'..','simulations',str(x),'results/thickness.csv')
			data, zd = read_vahid_file(fl)
			#data_pd=pd.read_csv(fl,names=variables)
			
			#print(data.thickness.iloc[-1])
			#print(data.tip_vel.iloc[-1])

			features1.append(data.thickness.iloc[-1])
			features2.append(data.tip_vel.iloc[-1])
			features3.append(data.TIME.iloc[-1])

			combined_last_rows=combined_last_rows.append(data.iloc[-1,:])
		except Exception as error:
			print("An exception occurred:", error) # An exception occurred: division by zero
			#print(lg)
			#log = log.append(lg)

        # Update the progress bar
		progress_bar.update(1)        

	# Close the progress bar
	progress_bar.close() 

	
	print('')
	print(' ******** example loaded result file ******** ')
	print(data)

	data = {'features1': features1,
			'features2': features2,
			'features3': features3
			}
	## Create PD dataframe frpm features
	features_df = pd.DataFrame(data)

	print('')
	print(' ******** Extracted features **************** ')
	print(features_df)

	## 
	print('')
	print(' ******** inputs and outputs are combined ******** ')
	input_outputs = pd.concat([samples,features_df],axis=1)


	input_outputs['features1'] = input_outputs['features1'] - input_outputs['features1'].iloc[0]
	data = input_outputs.dropna()

	#print(log)
	features_df.to_csv('../data/features_df.csv',na_rep="",index=False)
	combined_last_rows.to_csv('../data/combined_last_rows.csv',na_rep="",index=False)
	input_outputs.to_csv('../data/input_outputs_df.csv',na_rep="",index=False)
	
	
	plot_data(data)
	
	
