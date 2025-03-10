#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 18:42:54 2023

@author: attari.v
"""

import logging
#from pyDOE2 import *
from SALib.sample import latin
import pandas as pd
import numpy as np

##

from SALib.sample import latin

def _lhs_init(par_names, bounds, samples):
    """
    Returns LHS samples.

    :param par_names: List of parameter names
    :type par_names: list(str)
    :param bounds: List of lower/upper bounds,
            must be of the same length as par_names
    :type bounds: list(tuple(float, float))
    :param int samples: Number of samples
    :return: DataFrame
    """
    
    print('\n')
    print(' **************** Generating samples - LHS ***************** \n')
    
    # Define the problem
    problem = {
        'num_vars': len(par_names),
        'names': par_names,
        'bounds': bounds
    }
        
    # Generate LHS samples
    lhs_samps = latin.sample(problem, samples)
        
    # Convert samples to a dictionary
    par_vals = {}
    for i, par in enumerate(par_names):
        par_vals[par] = lhs_samps[:, i]
        
    # Convert dict(str: np.ndarray) to pd.DataFrame
    par_df = pd.DataFrame(columns=par_names, index=np.arange(samples))
    for i in range(samples):
        for p in par_names:
            par_df.loc[i, p] = par_vals[p][i]
            
    print(par_df)
            
    return par_df

#def _lhs_init(par_names, bounds, samples, criterion='c'):
#   """
#   Returns LHS samples.
#
#   :param par_names: List of parameter names
#   :type par_names: list(str)
#   :param bounds: List of lower/upper bounds,
#           must be of the same length as par_names
#   :type bounds: list(tuple(float, float))
#   :param int samples: Number of samples
#   :param str criterion: A string that tells lhs how to sample the
#           points. See docs for pyDOE.lhs().
#   :return: DataFrame
#   """
#
#   print(' \n');
#   print(' **************** Generating samples - LHS ***************** \n');
#
#   lhs_samps = lhs(len(par_names), samples=samples, criterion='c')
#   par_vals = {}
#   for par, i in zip(par_names, range(len(par_names))):
#       par_min = bounds[i][0]
#       par_max = bounds[i][1]
#       par_vals[par] = lhs_samps[:, i] * (par_max - par_min) + par_min
#
#   # Convert dict(str: np.ndarray) to pd.DataFrame
#   par_df = pd.DataFrame(columns=par_names, index=np.arange(samples))
#   for i in range(samples):
#       for p in par_names:
#           par_df.loc[i, p] = par_vals[p][i]
#
#   #logger = logging.getLogger(GA.__name__)
#   #logger.info('Initial guess based on LHS:\n{}'.format(par_df))
#   return par_df


from SALib.sample import sobol_sequence

def _sobol_init(par_names, bounds, samples):
    """
    Returns Sobol samples.

    :param par_names: List of parameter names
    :type par_names: list(str)
    :param bounds: List of lower/upper bounds,
            must be of the same length as par_names
    :type bounds: list(tuple(float, float))
    :param int samples: Number of samples
    :return: DataFrame
    """
    
    print('\n')
    print(' **************** Generating samples - Sobol ***************** \n')
    
    num_vars = len(par_names)
    
    # Generate Sobol samples
    sobol_samps = sobol_sequence.sample(samples, num_vars)
    
    par_vals = {}
    for par, i in zip(par_names, range(len(par_names))):
        par_min = bounds[i][0]
        par_max = bounds[i][1]
        par_vals[par] = sobol_samps[:, i] * (par_max - par_min) + par_min
        
    # Convert dict(str: np.ndarray) to pd.DataFrame
    par_df = pd.DataFrame(columns=par_names, index=np.arange(samples))
    for i in range(samples):
        for p in par_names:
            par_df.loc[i, p] = par_vals[p][i]
            
    return par_df

    
# ## GSD design 

# levels = [3, 3, 3, 3]
# reduction = 2       # Reduce the number of experiment to approximately a third.

# # ## create the levels
# samps = _gsd_init(par_names, levels, reduction)
# print(len(samps))

# ## name the columns for record 

# df = pd.DataFrame(samps, columns=['Ls','Le','W','eps'])

# df['W'] = np.where(df['W'] == 0, 2e-1, df['W'])
# df['W'] = np.where(df['W'] == 1, 4e-1, df['W'])
# df['W'] = np.where(df['W'] == 2, 6e-1, df['W'])

# df['eps'] = np.where(df['eps'] == 0, 1e-2, df['eps'])
# df['eps'] = np.where(df['eps'] == 1, 5e-2, df['eps'])
# df['eps'] = np.where(df['eps'] == 2, 5e-3, df['eps'])

# df['alp'] = np.where(df['alp'] == 0, 0.45, df['alp'])
# df['alp'] = np.where(df['alp'] == 1, 0.50, df['alp'])
# df['alp'] = np.where(df['alp'] == 2, 0.55, df['alp'])

# #df['bet'] = np.where(df['bet'] == 0, 0.40, df['bet'])
# #df['bet'] = np.where(df['bet'] == 1, 0.50, df['bet'])
# #df['bet'] = np.where(df['bet'] == 2, 0.60, df['bet'])

# df.to_csv('samples_GSD3.csv',index=False)
# df