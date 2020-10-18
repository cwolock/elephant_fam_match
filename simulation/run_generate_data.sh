#!/bin/bash

# num_combos is the number of parameter combinations (see generate_data.R)
# 8 loci numbers * 2 species * 4 relationships = 56
# 1 is the sim name
# 2 is the number of reps per combo (1)
# 3 is the number of reps performed per job (1)

num_combos=56
njobs=`expr $2 / $3 \* $num_combos`

qsub -cwd -e iotrash/ -o iotrash/ -t 1-$njobs ./call_generate_data.sh $1 $2 $3
