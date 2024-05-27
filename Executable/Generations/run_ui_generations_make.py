#!/usr/bin/env python

from os import system

log = f'run_ui_generations.log'

# Explanations follow each option as a comment after '#', with (defaults) in parentheses. 
O = ' -o ../../Output/Generations/generations.csv' # -o output filename
I = ' -i ../../Data/Generations/generations0.csv' # -i input filename' -w ../../Data/3_add_watchers_1.csv'
C = ' -c 20' # the maximum number of generations in the cdf output.

system( f'python run_ui_generations.py {O} {I} {C} > {log}' )

