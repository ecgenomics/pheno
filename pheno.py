#!/usr/bin/env python3

#        _                      
#       | |                     
#  _ __ | |__   ___ _ __   ___  
# | '_ \| '_ \ / _ | '_ \ / _ \ 
# | |_) | | | |  __| | | | (_) |
# | .__/|_| |_|\___|_| |_|\___/ 
# | |                           
# |_|                           
#
# Phenotype simulation tool

# Core script.

from modules.functions import *

### Testing area

tree = Phylo.read("tree.nw", "newick")

simwrapper("bm", 10, "tree.nw", "BM")
simwrapper(sw_name = "ou", sw_number_of_simulations = 10, sw_input_tree = "tree.nw", sw_mode = "OU", sw_seed = 0, sw_theta = 0.5, sw_mu = 8, sw_sigma = 0.22)