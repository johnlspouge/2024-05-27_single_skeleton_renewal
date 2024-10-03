#!/usr/bin/env python

from os import system

log = 'run_fecundity_negative_binomial.log'

# Explanations follow each option as a comment after '#', with (defaults) in parentheses. 
O = ' -o ../../Output/Fecundity/' # -o output directory name
K = ' -k 0.1 10.0' # K_ITER gives the difference and upper bound for iteration over negative binomial dispersion k.
P = ' -p 1.1 0.001' # P_ITER gives the factor and lower bound for iteration over negative binomial success probability p.

system( f'python run_fecundity_negative_binomial.py {O} {K} {P} > {log}' )

