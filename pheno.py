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

simulate_trait(
        
        input_tree = "tree.nw",
        mode = "bm"
                    )

simulate_trait(
        
        input_tree = "tree.nw",
        mode = "ou",
        st_mu = 4,
        st_dt = 0.1,
        st_sigma = 0.002,
        st_theta = 0.02
                    )

#print(simBM(1,0.1,4))

### End testing area

