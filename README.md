# Elephant familial matching

This repo contains the code used for familial matching and related analyses in the manuscript "Familial matching of tusks delineates the size and connectivity of transnational criminal organizations."

The code is organized into two main sections.

## Simulation code 

This code was used to perform the simulations found in the paper. Using elephant reference data, we artifically create relatives of various degrees (parent-offspring, full sibs, etc.) as well as unrelated individuals. We generate familial matching likelihood ratio distributions for true relatives and for unrelated individuals, and compare these distributions in order to choose a likelihood ratio cutoff for determining a match. 

This code is meant to be run in a cluster environment. 

## Data analysis code

This code was used to perform familial matching on the tusks analyzed in the paper. For each pair of tusks, we calculate the likelihood ratios for all possible relationships. There are also useful scripts for post-processing, including applying the simulaton-determined likelihood ratio threshold, and analyzing the number of matches between seizures.

The script `calculate_LRs.R` calculates pairwise LRs between all samples in the new sample file provided, plus LRs between pairs of samples where one is new and the other old. It does not calculate LRs between two old samples. You should run this script from inside the directory where it lives, as it relies on helper functions in `LR_functions.R`.  
