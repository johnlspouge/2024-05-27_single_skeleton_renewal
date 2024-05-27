#!/usr/bin/env python
"""
Calculates the statistics of single skeleton Galton-Watson process for gamma-distributed latent and infectious periods. 
"""
import sys
sys.path.insert(0,"../modules")

import argparse
from os.path import isfile, exists, dirname
from os import mkdir
from io import StringIO

import numpy as np
import pandas as pd
from math import log

from jls_branching_process import Branching_Process_Factory
from jls_epidemic_exponent import theta_solve

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
    cdf_max = argument.cdf_max
    df_cdf = pd.DataFrame(columns=[i for i in range(cdf_max+1)])
    lambdas = [] # exponential rate of infection
    doubling_times = []  # doubling time for infections
    ks = [] # dispersion in negative binomial offspring distribution
    ps = []  # probabiity in negative binomial offspring distribution
    qs = [] # extinction probability
    gammas = []  # renewal probability = mean offspring in doomed lineage
    geom_means = [] # mean total in doomed lineage
    geom_st_devs = []  # st dev total in doomed lineage
    # Iterates through rows on the input DataFrame.
    for index, row in df.iterrows(): 
        (e_mu, e_kappa, i_mu, i_kappa, r0) = row.to_list()
        print('parameter set', index, ':')
        print('    Latent Gamma(', e_mu, e_kappa, ')')
        print('    Infectious Gamma(', i_mu, i_kappa, ')')
        print('    R_0 (', r0, ')')
        # exponential rate of infection (lambda)
        theta = theta_solve(e_mu, e_kappa, i_mu, i_kappa, r0) # exponential rate (lambda)
        lambdas.append(theta)
        doubling_times.append(log(2.0)/theta)
        # branching process summary statistics
        bp = Branching_Process_Factory(r0=r0, dispersion=i_kappa) # Negative Binomial
        ks.append(i_kappa)
        ps.append(i_kappa/(i_kappa+r0)) 
        qs.append(bp.q())
        gamma = bp.gamma()
        gammas.append(gamma)
        geom_means.append(gamma/(1.0-gamma))
        geom_st_devs.append(gamma**0.5/(1.0-gamma))
        # skeleton branching process cdf G
        df_cdf.loc[index] = [1.0-gamma**(i+1) for i in range(cdf_max+1)]
    df = argument.df[COLS]
    df['lambda'] = np.array(lambdas) # exponential rate of infection
    df['doubling_time'] = np.array(doubling_times) # doubling time for infections
    df['k'] = np.array(ks) # dispersion in negative binomial offspring distribution
    df['p'] = np.array(ps)  # probabiity in negative binomial offspring distribution
    df['q'] = np.array(qs)  # probabiity in negative binomial offspring distribution
    qs = [] # extinction probability
    df['gamma'] = np.array(gammas)  # renewal probability = mean offspring in doomed lineage
    df['geom_mean'] = np.array(geom_means) # mean descendants in doomed lineage
    df['geom_st_dev'] = np.array(geom_st_devs) # st dev descendants in doomed lineage
    df = pd.concat([df,df_cdf], axis=1)
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
    if not isinstance(argument.cdf_max, int) or argument.cdf_max <= 0:
        raise ValueError(f'argument.cdf_max "{argument.cdf_max}" must be a positive integer.')
        
def getArguments():
    parser = argparse.ArgumentParser(description='Calculates the exponential growth lambda for SEIR model, where E and I are gamma-distributed.\n')
    parser.add_argument("-o", "--ofn", dest="ofn", type=str, default='../../Output/Exponent/exponent.csv', 
                        help="OFN contains the output DataFrame with the exponential growth lambda.", metavar="OFN")
    parser.add_argument("-i", "--ifn", dest="ifn", type=str, required=True,  
                        help="IFN is a CSV with DataFrame whose columns define the parameters of the gamma distributions.", metavar="IFN")
    parser.add_argument("-c", "--cdf_max", dest="cdf_max", type=int, required=True,  
                        help="CDF_MAX counts the atoms in the cdf of G, P(G<=g).", metavar="CDF_MAX")
    return parser
    
if __name__ == "__main__":
    main()
