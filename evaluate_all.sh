#!/usr/bin/env bash

# This bash script can be used to evaluate the generated
# modules for mexcowalk model as well as other models.
# The required python libraries are numpy, scipy, pandas,
# networkx and tqdm. The models to eveluate must be saved
# in a folder with the name of the model in the following
# order:
#########################################################
#
#		out/connected_components_isolarge/model_name
#
#########################################################


cd src


#########################################################
#
#	Evaluate Given Models
#
#########################################################


# following arguments can be passed onto this python file
#
# -m	str			Name of the model | default 'mexcowalk'
#					Add more models by adding space and
#					followed by model name i.e 'model1' 'model2' 'model3'
# -ns	int			Maximum number of genes | default 1000
# -ne	int			Minimum number of genes | default 100
# -step	int			Stepsize of decrement from num_start to
#					num_end | default 100
#
#########################################################


printf "\nRunning Evaluations...\n\n"

python evaluations.py \
	-m		'mexcowalk' \
	-ns 	1000 \
	-ne 	100 \
	-step	100

printf "\nDONE!!!\n"
