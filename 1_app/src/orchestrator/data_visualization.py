#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 08:47:18 2023

@author: attari.v
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def load_data(fl):
    
    df = pd.read_csv(fl)
    
    return df


fl = '../data/input_outputs_df.csv'
data = load_data(fl)

#data['features1'] #= data['features1'] #- .1874344E-03
#data['features1']
#data = data.dropna()
#data = data[data['features1'] > 0]
#data = data[data['features2'] > 0]
#data = data[data['features1'] > 1e-5]
#data = data[data['features1'] < 5e-5]

#data = data[data['features1'] < 2.0015e-5] # upper level 
#data = data[data['features1'] > 1.5e-5]

#data['features6'] = data['features4'] + data['features5']

#data_sorted = data.sort_values(by=['features4'])

print(data)
#print(data_sorted)

fl = '../data/combined_last_rows.csv'
combined_last_rows = load_data(fl)

fl = '../data/features_df.csv'
features_df = load_data(fl)

#fig.figure()
# p1 = pd.plotting.scatter_matrix(data, alpha=0.2, figsize=(15, 12), grid=True)
# for i, axs in enumerate(p1):
#     for j, ax in enumerate(axs):
#         if i != j:  # only the scatter plots
#             ax.set_xscale('log')
#             ax.set_yscale('log')


# # plot scatter matrix using seaborn
# sns.set_theme(style="ticks")
# pp = sns.pairplot(data=data)

# pp.fig.set_figheight(20)
# pp.fig.set_figwidth(20)

# log_xcolumns = data.columns.to_list()
# log_ycolumns = ["Ls","Le","W","eps","features1"]

# for ax in pp.axes.flat:
#     if ax.get_xlabel() in log_xcolumns:
#         ax.set(xscale="log")
#     if ax.get_ylabel() in log_ycolumns:
#         ax.set(yscale="log")        
        

########

# %%
        
log_xcolumns = ["Ls","Le","W","eps"]
log_ycolumns = ["features1"]

lables_xcolumns = ["Interface Mobility","Reaction Rate","Barrier Height of Transformation","Gradient Energy Coef."]
lables_ycolumns = ["features1: Thickness (m/s)"]


time = [60, 136, 215]; # min
time = [i * 60 for i in time]  # seconds
V = np.array([1.16e-6, 0.76e-6, 0.65e-6])*1e-2; #cm/s --> m/s

M = 24.305/1000   #g/mol  /1000 ---> Kg/mol
D = 1.03e-7/10000 #cm2/s  /10000 ---> m^2/s
c = 0.1/1000*1e6   #M -> /1000 ---> mol/cm^3 * 1e6 ---> mol/m^3
ro= 1.74*1000     #g/cm³  *1000 ---> Kg/m^3

t               = 60*60
theory_R_250min = ((2*M*D*c*t)/(ro))**0.5
g_truth_250min  = np.repeat(theory_R_250min,np.size(data,0))

t               = 30*60
theory_R_60min = ((2*M*D*c*t)/(ro))**0.5

g_truth_60min = np.repeat(theory_R_60min,np.size(data,0))


#i=1
j=0
fig, ax = plt.subplots(1,4,sharey=True, figsize = (24,6))
for idx_y in log_ycolumns:
    for idx_x in log_xcolumns:
        print(j)
        print(idx_x)
        scatter = ax[j].scatter(data[idx_x],data[idx_y], alpha=0.7, edgecolors="k", c=data.features3, s=data.features3/60)
        scatter2= ax[j].plot(data[idx_x],g_truth_250min, "--", c='black', label='250 min')
        scatter3= ax[j].plot(data[idx_x],g_truth_60min , "-", c='black', label='60 min')
        ax[j].set_xscale("log");
        ax[j].set_yscale("log");
        ax[j].set_xlabel(lables_xcolumns[j])
        ax[j].set_ylabel(idx_y+' :Electrode Thickness (m)')
        #ax.fontsize = 15

        # produce a legend with a cross-section of sizes from the scatter
        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
        legend2 = ax[j].legend(handles, labels, loc="upper right", title="Time (min)")        
        
        j = j+1

fig.savefig('data_feature1.jpg',dpi=200)

log_xcolumns = ["Ls","Le","W","eps"]
log_ycolumns = ["features2"]


pause
# %%

j=0
fig, ax = plt.subplots(1,4,sharey=True, figsize = (24,6))
for idx_y in log_ycolumns:
    for idx_x in log_xcolumns:
        print(idx_x)
        scatter = ax[j].scatter(data[idx_x],data[idx_y], alpha=0.7, edgecolors="k", c=data.features3, s=data.features3/60)
        ax[j].set_xscale("log");
        ax[j].set_yscale("log");
        ax[j].set_xlabel(lables_xcolumns[j])
        ax[j].set_ylabel(idx_y+str(' :Interface Velocity (m/s)'))
        #ax.fontsize = 15

        # produce a legend with a cross-section of sizes from the scatter
        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
        legend2 = ax[j].legend(handles, labels, loc="upper right", title="Time (min)")        
        
        j = j+1
        
fig.savefig('data_feature2.jpg',dpi=200)

log_xcolumns = ["Ls","Le","W","eps"]
log_ycolumns = ["features3"]

# %%

j=0
fig, ax = plt.subplots(1,4,sharey=True, figsize = (24,6))
for idx_y in log_ycolumns:
    for idx_x in log_xcolumns:
        print(idx_x)
        scatter = ax[j].scatter(data[idx_x],data[idx_y], alpha=0.7, edgecolors="k", c=data.features3, s=data.features3/60)
        ax[j].set_xscale("log");
        ax[j].set_yscale("log");
        ax[j].set_xlabel(lables_xcolumns[j])
        ax[j].set_ylabel(idx_y+str(' :Time (s)'))
        #ax.fontsize = 15
        
        # produce a legend with a cross-section of sizes from the scatter
        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.4)
        legend2 = ax[j].legend(handles, labels, loc="upper right", title="Time (min)")        

        j = j+1
        
fig.savefig('data_feature3.jpg',dpi=200)

# %%


## Find useful data:
## data_useful = data.loc[(data['features2'] < 1e-7 ) & (data['features2'] != 0) & (data['features2'] > 0) ]

log_xcolumns = ["Ls","Le","W","eps"]
log_ycolumns = ["Ls","Le","W","eps"]

i=0
fig, ax = plt.subplots(4,4, figsize = (20,20))
for lb_x, idx_x in zip(lables_xcolumns,log_xcolumns):
    j=0
    for lb_y, idx_y in zip(lables_xcolumns,log_xcolumns):
        print(idx_x)
        ax[i,j].scatter(data[idx_x],data[idx_y], alpha=0.7, edgecolors="k", c=data.features1)
        ax[i,j].set_xscale("log");
        ax[i,j].set_yscale("log");
        ax[i,j].set_xlabel(lb_x)
        ax[i,j].set_ylabel(lb_y)
        #ax.fontsize = 15
        
        j = j+1
    i = i+1
        
#fig.savefig('data_feature2.jpg',dpi=200)

# %%

##### Python code to reproduce Fig. 4(a):

fig, ax = plt.subplots( figsize = (8,8))

exp_time = [60, 136, 215]; # min
exp_time = np.array(exp_time)*60 # seconds
exp_V = np.array([1.16e-6, 0.76e-6, 0.65e-6])*1e-2; #cm/s *1e-2  ---> m/s

M = 24.305/1000   #g/mol  /1000 ---> Kg/mol
D = 1.03e-5/10000 #cm2/s  /10000 ---> m^2/s
c = 0.1/1000*1e6   #M -> /1000 ---> mol/cm^3 * 1e6 ---> mol/m^3
ro= 1.74*1000     #g/cm³  *1000 ---> Kg/m^3

t        = np.linspace(4,250,100)*60
theory_V = ((M*D*c)/(2*ro))**0.5*(1/(t)**0.5) # m/s

ax.plot(t,theory_V,'--')
ax.plot(exp_time, exp_V, 's')

ax.set_xlabel('Time [seconds]')
ax.set_ylabel('Theoretical rate [m/s]')

# %%


## plot 1
## plot 1
## plot 1 ------> velocity, m/s
## plot 1
## plot 1

M = 24.305 #g/mol
D = 1.03e-7 #1.03e-7 #cm2/s
c = 0.1/1000 #M -> mol/cm^3
ro= 1.74   #g/cm³

t = np.linspace(4,220,100) # min
theory_V = ((M*D*c)/(2*ro))**0.5*(1/(t)**0.5)*1e-2 # cm/s--> m/s

fig, ax = plt.subplots( figsize = (8,8))

time = [64, 136, 215]; # min
V = np.array([1.16e-6, 0.76e-6, 0.65e-6])*1e-2; #cm/s --> m/s

t_1hr = 60 # min
theory_R_1hr = ((2*M*D*c*t_1hr)/(ro))**0.5*1e-2

ax.plot(t,theory_V,'--')
#ax.plot(t_1hr,theory_R_1hr,'-s',c='black')
#ax.text(t[-1], theory_V[-1], theory_V[-1], fontsize=12)
ax.plot(time, V, 's')

ax.set_xlabel('Time [min]')
ax.set_ylabel('Theoretical velocity [m/s]')

## plot 1
## plot 1
## plot 1 ------> velocity, paper cm/s
## plot 1
## plot 1

M = 24.305 #g/mol
D = 1.03e-7 #1.03e-7 #cm2/s
c = 0.1/1000 #M -> mol/cm^3
ro= 1.74   #g/cm³

t = np.linspace(4,220,100) # min
theory_V = ((M*D*c)/(2*ro))**0.5*(1/(t)**0.5)

fig, ax = plt.subplots( figsize = (8,8))

time = [64, 136, 215]; # min
V = np.array([1.16e-6, 0.76e-6, 0.65e-6])*1e6; #cm/s

t_1hr = 60 # min
theory_R_1hr = ((2*M*D*c*t_1hr)/(ro))**0.5

ax.plot(t,theory_V*1e6,'--')
#ax.plot(t_1hr,theory_R_1hr,'-s',c='black')
#ax.text(t[-1], theory_V[-1], theory_V[-1], fontsize=12)
ax.plot(time, V, 's')

ax.set_xlabel('Time [min]')
ax.set_ylabel('Theoretical [x10$^{-6}$ cm/s]')

## Anson equation

slope = 1.7e-6 # C/s^1/2
r = 0.4 # 4 mm to cm
F = 96485 # C/mol
c = 0.1/1000 # mol/cm^3
D = (slope/(2*F*3.14*r**2*c*3.14**(-0.5) ))**2
print('D [cm^2/s]=',D)

# %% 

## plot 1
## plot 1
## plot 1 ------> thickness
## plot 1
## plot 1

fig, ax = plt.subplots( figsize = (8,8))

time = [i * 60 for i in time]  # seconds

time = np.array([60, 136, 215])*60; # min
V = np.array([1.16e-6*60*60, 0.76e-6*136*60, 0.65e-6*215*60])*1e-2*1e6; #cm/s

M = 24.305/1000   #g/mol  /1000 ---> Kg/mol
D = 1.03e-7/10000 #cm2/s  /10000 ---> m^2/s
#D = 1.03e-7/10000 #cm2/s  /10000 ---> m^2/s
c = 0.1/1000*1e6   #M -> /1000 ---> mol/cm^3 * 1e6 ---> mol/m^3
ro= 1.74*1000     #g/cm³  *1000 ---> Kg/m^3

t        = np.linspace(4,250,100)*60
theory_R = ((2*M*D*c*t)/(ro))**0.5*1e6

t_1hr = 60*60 # [s]
theory_R_1hr = ((2*M*D*c*t_1hr)/(ro))**0.5*1e6 # 1e6 convert m to um

ax.plot(t,theory_R,'--')
ax.plot(t_1hr,theory_R_1hr,'-s',c='black')
#ax.text(t[-1], theory_V[-1], theory_V[-1], fontsize=12)
ax.text(t_1hr, theory_R_1hr, theory_R_1hr, fontsize=12)
ax.plot(time, V, 's')

ax.set_xlabel('Time [s]')
ax.set_ylabel('Theoretical Radius [$\mu$ m]')
ax.set_yscale('log')

print(theory_R[-1])

# %%
## Mg self-diffusion

fig, ax = plt.subplots( figsize = (8,8))

R = 8.314
T = np.linspace(300,600,20) # K
D_Mg_hcp = 2.9e-5*np.exp(-125748/(R*T))

ax.plot(T,D_Mg_hcp)
ax.text(T[0], D_Mg_hcp[0], D_Mg_hcp[0], fontsize=12)

ax.set_xlabel('Temperature [K]')
ax.set_ylabel('D [m2/s]')
ax.set_yscale('log')

M = D_Mg_hcp[0]/(R*T[0])

print('Mobility:=',M)

