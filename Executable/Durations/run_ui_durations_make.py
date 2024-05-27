#!/usr/bin/env python

from os import system

log = f'run_ui_durations.log'

# Explanations follow each option as a comment after '#', with (defaults) in parentheses. 
O = ' -o ../../Output/Durations/durations.csv' # -o output filename
I = ' -i ../../Data/Durations/durations0.csv' # -i input filename' -w ../../Data/3_add_watchers_1.csv'
R = ' -r 1000' # the realizations of single skeleton renewal simulation

system( f'python run_ui_durations.py {O} {I} {R} > {log}' )

