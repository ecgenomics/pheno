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

# Functions script

import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

### FUNCTION simOU()
### Simulates an Ornstein-Uhlenbeck process, starting from a
### known seed value (s = 0 by default)

def simOU(theta, mu, sigma, t, dt, s = 0):

    num_steps = int(t/dt)

    # Init result vector x
    x = np.zeros(num_steps)
    x[0] = s

    # White noise generation
    dW = np.random.normal(0, np.sqrt(dt), size=num_steps)

    # Process simulation
    for t in range(1, num_steps):
        x[t] = x[t - 1] + theta * (mu - x[t - 1]) * dt + sigma * dW[t]
    
    return x

### FUNCTION simBM()
### Simulates a Bwonian Motion process, starting from a
### known seed value (s = 0 by default)

def simBM(t, dt, s=0):

    num_steps = int(t/dt)

    # White noise generation
    dW = np.random.normal(0, np.sqrt(dt), size=num_steps)
    dW = np.insert(dW,0, s, axis =0)
    
    # Wiener process simulation
    W = np.cumsum(dW)
    
    return W

### FUNCTION cdist()
### Calculates the distance between two nodes

def cdist(tree, node1, node2):
    # Trova i percorsi dai due nodi all'antenato comune più vicino
    path1 = tree.trace(node1, node2)
    path2 = tree.trace(node2, node1)
    
    # Calcola la distanza sommando le lunghezze dei rami dei percorsi
    distance = sum(branch.length for branch in path1) + sum(branch.length for branch in path2)
    
    return distance

### FUNCTION simulate_trait()
### simulates trait evolution on a tree

def simulate_trait(
        
        input_tree,
        mode,
        seed = 0,
        st_theta = 0,
        st_mu = 0,
        st_sigma = 0,
        st_dt = 0.1,
        tree_format = "newick"

                    ):
    
    class node_values():

        def __init__(self):
            self.name = ""
            self.sim_val = 0

            self.children = []
            self.terminals = []
            self.leaves = []
            self.path2child = {}
            self.child_val = {}

    tree = Phylo.read(input_tree,tree_format)
    
    # Initialize the tree (give a name to internal nodes int_#node)

    counter = 1

    for n in tree.find_clades():
        if n.name is None:
            n.name = f"int_{counter}"
            counter += 1
    
    # Collect data
    
    child_dict = {}
    path_dict = {}
    vals_dict = {}
    
    # Iterate the simulation
    for clade in tree.find_clades(order='level'):
        z = node_values()
        z.name = clade.name

        try:
            z.sim_val = z.vals_dict[z.name]
        except:
            
            z.sim_val = seed

        if clade.clades:
            z.children = [child.name for child in clade.clades]
            child_dict[z.name] = [child.name for child in clade.clades]

            z.terminals = [t.name for t in clade.get_terminals()]
        
        print(z.name, z.sim_val, z.children, z.terminals)

        if clade.clades:
            print(clade.distance(clade.clades[0]))


    Phylo.draw(tree)