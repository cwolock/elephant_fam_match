#!/bin/bash

Rscript ../../scripts/cluster/generate_data.R --sim-name $1 --nreps-total $2 --nreps-per-job $3
