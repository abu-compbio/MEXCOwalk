#!/usr/bin/env bash

# This bash script can be used to run all python files
# sequentially for the mexcowalk algorithm. The required 
# python libraries are numpy, scipy and networkx. The 
# input files are provided within the data folder. If 
# the user wishes to use a different input file, the 
# filenames in the utils.py file need to be changes.

cd src


#########################################################
#
#	Compute Edge Weights
#
#########################################################

# following arguments can be passed onto this python file
#
# -m	str			Name of the model | default 'mexcowalk'
# -t	float	 	Mutex threshold | default 0.7
# -v				Verbose, prints additional data during 
#					runtime
#
#########################################################

printf "Computing Edge Weights...\n\n"

python compute_edge_weights.py \
	-v


#########################################################
#
#	Perform Random Walk
#
#########################################################

# following arguments can be passed onto this python file
#
# -m	str			Name of the model | default 'mexcowalk'
# -s			 	Stores all random walk matrices when 
#					called
# -v				Verbose, prints additional data during 
#					runtime
#
#########################################################

printf "\nRunning Random Walk...\n\n"

python random_walk.py \
	-s	\
	-v			


#########################################################
#
#	Find Connected Components and Split
#
#########################################################

# following arguments can be passed onto this python file
#
# -m	str			Name of the model | default 'mexcowalk'
# -k	int			Minimum length of a module | default 3
# -ns	int			Maximum number of genes | default 1000
# -ne	int			Minimum number of genes | default 100
# -step	int			Stepsize of decrement from num_start to
#					num_end | default 100
# -ts	float		Starting threshold | default 0.0002 for
#					1000 genes
# -cc			 	Outputs modules before splitting when 
#					called
# -v				Verbose, prints additional data during 
#					runtime
#
#########################################################


# Change -ts to 0.00008 and ns to 2500 in order to find 
# top 2500 genes. -ts to 0.0002 and ns to 1000 for top 
# 1000 genes.
#
#########################################################

printf "\nFinding Connected Components...\n\n"

python connected_components_isolarge.py \
	-k		3 \
	-ns 	1000 \
	-ne 	100 \
	-step	100 \
	-ts		0.0002 \
	-v			

printf "\nDONE!!!\n"

