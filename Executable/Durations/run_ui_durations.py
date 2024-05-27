#!/usr/bin/env python
"""
Calculates the statistics of the duration of the single skeleton process for gamma-distributed latent and infectious periods. 
"""
import sys
sys.path.insert(0,"../modules")

import argparse
from os.path import isfile, exists, dirname
from os import mkdir
from io import StringIO

import numpy as np
import pandas as pd
from scipy.stats import geom, gamma, uniform

from jls_branching_process import Branching_Process_Factory

# names of the relevant columns
COLS = ['e->i_mean',
        'e->i_dispersion',
        'i->r_mean',
        'i->r_dispersion',
        's->e:i_R_0']

def main():
    parser = getArguments()
    argument = parser.parse_args()
    check(argument) 
    df = argument.df
    realization_num = argument.realization_num
    df_realization = pd.DataFrame(columns=[i for i in range(realization_num+1)])
    # Iterates through rows on the input DataFrame.
    for index, row in df.iterrows(): 
        (e_mu, e_kappa, i_mu, i_kappa, r0) = row.to_list()
        print('parameter set', index, ':')
        print('    Latent Gamma(', e_mu, e_kappa, ')')
        print('    Infectious Gamma(', i_mu, i_kappa, ')')
        print('    R_0 (', r0, ')')
        durations = []
        # branching process summary statistics
        bp = Branching_Process_Factory(r0=r0, dispersion=i_kappa) # Negative Binomial
        q = bp.q()
        gamma_bp = bp.gamma()
        for r in range(realization_num):
            # duration of single skeleton renewal
            g = geom.rvs(1.0-gamma_bp, size=1)[0]
            if g == 1:
                durations.append(0.0)
                continue
            latent = gamma.rvs((g-1)*e_kappa, scale=e_mu/e_kappa, size=1)
            infects = gamma.rvs(i_kappa+1.0, scale=i_mu/((1.0-q)*r0+i_kappa), size=g-1)
            uniforms = uniform.rvs(size=g-1)
            maternal_birth = 0.0
            for i in range(g-1):
                maternal_birth += infects[i]*uniforms[i]
            duration = latent+maternal_birth
            duration = duration.astype(float)[0]
            durations.append(duration)
        durations = [0.0]+sorted(durations)
        print('     sample_mean =', sum(durations)/realization_num)
        mean_maternal_birth = 0.5*i_mu*(i_kappa+1.0)/((1.0-q)*r0+i_kappa)
        print('            mean =', (e_mu+mean_maternal_birth)*gamma_bp/(1.0-gamma_bp))
        df_realization.loc[index] = durations
    df = argument.df[COLS]
    df = pd.concat([df,df_realization], axis=1)
    df.to_csv(argument.ofn, index=False)
# Reads the string from argument.ifn, a CSV DataFrame with columns COLS.
def to_df(string):
    ifh = StringIO(string)
    df = pd.read_csv(ifh, sep=',', dtype={COLS[0]:float,COLS[1]:float,COLS[2]:float,COLS[3]:float,COLS[4]:float})
    arr = df.to_numpy()
    if (np.abs(arr) <= 0.0).any():
        raise ValueError('The dataframe elements must be positive.')
    return df
# Checks arguments.    
def check(argument): 
    # argument.odir is the output directory.
    odir = dirname(argument.ofn)
    if not exists(odir):
        mkdir(odir)
    # argument.ifn contains a DataFrames with columns 
    #    'e->i_mean' : E gamma distribution mean
    #    'e->i_dispersion' : E gamma distribution dispersion
    #    'i->r_mean' : I gamma distribution mean
    #    'i->r_dispersion' : I gamma distribution mean
    #    's->e:i_R_0' : I mean basic reproduction number
    if not isfile(argument.ifn):
        raise ValueError(f'Input file "{argument.ifn}" does not exist.')
    try:
        with open(argument.ifn, 'r') as ifh:
            string = ifh.read()
    except:
        raise ValueError(f'The read of the input file "{argument.ifn}" failed.')
    try:
        argument.df = to_df(string)
    except:
        raise ValueError(f'The input file "{argument.ifn}" contained bad values.')
    if not isinstance(argument.realization_num, int) or argument.realization_num <= 0:
        raise ValueError(f'argument.realization_num "{argument.realization_num}" must be a positive integer.')
        
def getArguments():
    parser = argparse.ArgumentParser(description='Calculates the exponential growth lambda for SEIR model, where E and I are gamma-distributed.\n')
    parser.add_argument("-o", "--ofn", dest="ofn", type=str, default='../../Output/Exponent/exponent.csv', 
                        help="OFN contains the output DataFrame with the exponential growth lambda.", metavar="OFN")
    parser.add_argument("-i", "--ifn", dest="ifn", type=str, required=True,  
                        help="IFN is a CSV with DataFrame whose columns define the parameters of the gamma distributions.", metavar="IFN")
    parser.add_argument("-r", "--realization_num", dest="realization_num", type=int, required=True,  
                        help="REALIZATION_NUM counts the realizations of single skeleton renewal simulation.", metavar="REALIZATION_NUM")
    return parser
    
if __name__ == "__main__":
    main()
