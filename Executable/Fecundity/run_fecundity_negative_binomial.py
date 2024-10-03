#!/usr/bin/env python
"""
Calculates gamma, the mean offspring number of a dual subcritical GW process.
"""
import sys
sys.path.insert(0,"../modules")

import argparse
from os.path import exists
from os import mkdir

import numpy as np
import pandas as pd

from jls_branching_process import Negative_Binomial

def main():
    parser = getArguments()
    argument = parser.parse_args()
    check(argument) 
    k_difference, k_end = argument.k_iter # Increases k by k_difference from 0.0 until k_end < k.
    print('k iteration parameters :',argument.k_iter)
    p_factor, p_end = argument.p_iter # Decreases p by a factor of 1.0/p_factor from 1.0 until p < p_end.   
    print('p iteration parameters :',argument.p_iter)
    # Calculates df_q['k'] == df_gamma['k']. 
    ks = [] # dispersion of negative binomial
    k = 0.0
    while True:
        k += k_difference
        if k_end < k:
            break
        ks.append(k)
    # Initializes the columns of df_q & df_gamma whose heading is p.
    df_q = pd.DataFrame()
    df_q['k'] = np.array(ks) # mean offspring number
    df_gamma = pd.DataFrame()
    df_gamma['k'] = np.array(ks) # mean offspring number
    df_derivative = pd.DataFrame()
    df_derivative['k'] = np.array(ks) # mean offspring number
    # Calculates the columns of df whose heading is p.
    p = p_old = 1.0
    while True: 
        p /= p_factor
        if p < p_end:
            break
        qs = [] # renewal probability = extinction probability
        gammas = [] # renewal probability = mean offspring in doomed lineage
        derivatives = [] #
        k_old = 0.0
        q_old = gamma_old = 1.0
        for i,k in enumerate(ks):
            if k*(1.0-p)/p <= 1.0: # subcriticality
                qs.append(1.0)
                gammas.append(1.0)
                derivatives.append(0.0) # derivative of gamma wrt x
                continue
            assert isinstance(k, float) and 0.0 < k <= k_end
            assert isinstance(p, float) and p_end <= p < 1.0
            bp = Negative_Binomial(k, p) # Negative Binomial
            q = bp.q()
            assert isinstance(q, float)
            gamma = bp.gamma()
            assert isinstance(gamma, float)
            # Checks that values of q and gamma decrease.
            if p_old < 1.0 and q > df_q.at[i, p_old]:
                print('q p :', k, p, q, ';', df_q.at[i, p_old], k, p_old)
            if q > q_old:
                print('q k :', k, p, gamma, ';', gamma_old, k_old, p)
            if p_old < 1.0 and gamma > df_gamma.at[i, p_old]:
                print('gamma p :', k, p, gamma, ';', df_gamma.at[i, p_old], k, p_old)
            if gamma > gamma_old:
                print('gamma k :', k, p, gamma, ';', gamma_old, k_old, p)
            qs.append(q)
            gammas.append(gamma)
            # derivatives.append(-derivative(k, p, q, gamma)) # derivative of gamma wrt x
            k_old = k
            q_old = q
            gamma_old = gamma
        df_q[p] = np.array(qs) # extinction probability
        df_gamma[p] = np.array(gammas) # dual mean offspring number
        # df_derivative[p] = np.array(derivatives) # derivative of gamma wrt x
        p_old = p
    df_q.apply(pd.to_numeric, errors='raise')
    df_q.to_csv(argument.odir+'negative_binomial_q.csv', index=False)
    df_gamma.apply(pd.to_numeric, errors='raise')
    df_gamma.to_csv(argument.odir+'negative_binomial_gamma.csv', index=False)
    # df_derivative.apply(pd.to_numeric, errors='raise')
    # df_derivative.to_csv(argument.odir+'negative_binomial_negative_derivative.csv', index=False)
def derivative(k, p, q, gamma):
#    return 1.0-(1.0-p)/(1.0-gamma)*k*(1.0/p+(1.0-p)/(1.0-(1.0-p)*q)) # always positive
#    return (1.0-gamma)-(1.0-p)*((1.0-gamma)*q+k*(1.0+q-p)) # error
     return p*(1.0-(1.0-p)*q)*(1.0-gamma)-(1.0-p)*k*(1.0-(1.0-p)*q+p*(1.0-p)) # always positive
# Checks arguments.    
def check(argument): 
    # argument.odir is the output directory.
    if not argument.odir.endswith('/'):
        argument.odir += '/'
    odir = argument.odir
    if not exists(odir):
        mkdir(odir)
    # Checks argument.k_iter.
    if not isinstance(argument.k_iter, list) or len(argument.k_iter) != 2:
        raise ValueError(f'argument.k_iter "{argument.k_iter}" must be a list of length 2.')
    for i in range(2):
        x = argument.k_iter[i]
        if not isinstance(x, (int, float)) or x <= 0.0:
            raise ValueError(f'argument.k_iter "{x}" must be a positive int or a float.')
        argument.k_iter[i] = float(argument.k_iter[i])
    # Checks argument.p_iter.
    if not isinstance(argument.p_iter, list) or len(argument.p_iter) != 2:
        raise ValueError(f'argument.k_iter "{argument.p_iter}" must be a list of length 2.')
    for i in range(2):
        x = argument.p_iter[i]
        if not isinstance(x, (int, float)):
            raise ValueError(f'argument.p_iter "{x}" must be a float.')
        if i == 1 and (x <= 0.0 or 1.0 <= x):
            raise ValueError(f'argument.p_iter "{x}" must be a probability 0 < x < 1.')
        argument.p_iter[i] = float(argument.p_iter[i])
          
def getArguments():
    parser = argparse.ArgumentParser(description='Calculates gamma, the mean offspring number of a dual subcritical GW process.\n')
    parser.add_argument("-o", "--odir", dest="odir", type=str, default='../../Output/Fecundity/', 
                        help="ODIR is the output directory containing the output DataFrame-s with q and gamma, extinction probability and the dual mean offspring number.", metavar="OFN")
    parser.add_argument("-k", "--k_iter", dest="k_iter", nargs=2, type=float, required=True,  
                        help="K_ITER gives the difference and upper bound for k for iteration over negative binomial.", metavar="K_ITER")
    parser.add_argument("-p", "--p_iter", dest="p_iter", nargs=2, type=float, required=True,  
                        help="P_ITER gives the factor and lower bound for p for iteration over negative binomial p-s.", metavar="P_ITER")
    return parser
    
if __name__ == "__main__":
    main()
