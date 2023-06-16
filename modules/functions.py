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

import os
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

def simBM(t, dt, s):

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

### FUNCTION printdict()

def printdict(dictionary, out, tag):
    for x in dictionary.keys():
        print(tag, x, dictionary[x], file = out)

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
    
    ### CLASS node_values()
    ### Recaps all the simulation

    class simulation():

        def __init__(self):
            self.tag = "sim_"
            self.clades = []
            self.n2n_paths = {}
            self.node_values = {}
            self.leaf_values = {}


    ### CLASS node_values()
    ### Adds simulation features to the Clade() class
    ### Not derived from Clade(), but "parallel"

    class node_values():

        def __init__(self):
            self.name = ""
            self.sim_val = 0

            self.children = []
            self.terminals = []
            self.terminal_values = []
            self.terminal_values_mean = []
            self.terminal_values_sd = []

    # Read the tree
    tree = Phylo.read(input_tree,tree_format)
    

    # Create the simulation() object

    sim = simulation()

    # Initialize the tree (give a name to internal nodes int_#node)

    counter = 1

    for n in tree.find_clades():
        if n.name is None:
            n.name = f"int_{counter}"
            counter += 1
    
    # Iterate the simulation


    for clade in tree.find_clades(order='level'):
        z = node_values()
        z.name = clade.name
        try:
            z.sim_val = sim.node_values[clade.name]
        except:
            z.sim_val = seed
            sim.node_values[z.name] = z.sim_val
        
        if clade.clades:
            z.children = [child for child in clade.clades]

            for c in z.children:
                c_t = clade.distance(c)
                c_dt = st_dt
                seed = z.sim_val

                #### MODEL CHOICE ####

                if mode.lower() == "bm":

                    # BROWNIAN MOTION

                    out = simBM(c_t,c_dt,seed)

                elif mode.lower() == "ou":

                    # ORNSTEIN-ULLEHNBECK

                    out = simOU(mu = st_mu, sigma = st_sigma, theta = st_theta, t = c_t, dt = c_dt, s = seed)
                
                else:
                    print("incorrect model")
                    exit()

                last_value = list(out)[-1]
                sim.node_values[c.name] = last_value
                sim.n2n_paths[c.name] = out

            z.terminals = [t.name for t in clade.get_terminals()]
        else:
            sim.leaf_values[z.name] = z.sim_val

        sim.clades.append(z)
    
    # Work on the terminals

    for x in sim.clades:
        if len(x.terminals) > 0: 
            x.terminal_values = [sim.node_values[x] for x in x.terminals]
            x.terminal_values_mean = np.mean(x.terminal_values)
            x.terminal_values_sd = np.std(x.terminal_values)
    
    return sim

### FUNCTION simwrapper()
### repeats the simulation and saves the results in different files

def simwrapper(
        
        sw_name,                            # wrapper variables
        sw_number_of_simulations,

        sw_input_tree,                      # simulate_trait variables
        sw_mode,
        sw_seed = 0,
        sw_theta = 0,
        sw_mu = 0,
        sw_sigma = 0,
        sw_dt = 0.1,
        sw_tree_format = "newick",

        sw_rewrite = True                   # default wrapper variables

    ):    

    # Initialize the results folder

    results_folder = sw_name + "_folder"    

    try:
        os.mkdir(results_folder)
    except:
        if sw_rewrite == True:
            os.system("rm -r " + results_folder)
            os.mkdir(results_folder)
        else:
            pass
    # Counter
    simulation_counter = 0
    # Outputs
    nodevalues = open(results_folder + "/node.values", "w")
    leafvalues = open(results_folder + "/leaf.values", "w")
    paths = open(results_folder + "/paths.pernode", "w")
    clusters = open(results_folder + "/clusters.values", "w")



    while simulation_counter < sw_number_of_simulations:
        s = simulate_trait(     
        input_tree = sw_input_tree ,
        mode = sw_mode,
        seed = sw_seed,
        st_theta = sw_theta,
        st_mu = sw_mu,
        st_sigma = sw_sigma,
        st_dt = sw_dt,
        tree_format = sw_tree_format)

        # Counter
        simulation_counter += 1
        simtag = "cycle_" + str(simulation_counter)

        # Cluster file
        print("simulation", "node", "#leaves", "mean", "stdev", "leaf values", "mode= "+sw_mode, file=clusters)
        for x in s.clades:
            if len(x.terminal_values) > 0:
                print(simtag, x.name, len(x.terminal_values), x.terminal_values_mean, x.terminal_values_sd,x.terminal_values, file = clusters )
        print("#", file=clusters)

        # Paths file
        printdict(s.n2n_paths,paths, simtag)     
        print("# Mode =" + sw_mode, file=paths)
   
        printdict(s.leaf_values,leafvalues, simtag)  
        print("# Mode =" + sw_mode, file=leafvalues)
      
        printdict(s.node_values,nodevalues, simtag)        
        print("# Mode =" + sw_mode, file=nodevalues)

    # Output Closure
    clusters.close()
    paths.close()
    leafvalues.close()
    nodevalues.close()
    return 0