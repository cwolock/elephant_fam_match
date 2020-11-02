#!/usr/bin/env Rscript

#######################
### calculate_LRs.R ###
#######################
#
# calculates all pairwise familial matching LRs
# 
# command line arguments:
# 1) species: 'forest' or 'savannah'
# 2) reference file (.csv)  in long format, 2 lines per individual
#    rec is to include only individuals genotyped at >= 10 loci
# 3) sample file (.csv) in wide format, 1 line per individual
#
# output: tab-delimited file with pairwise LRs

# read in command line args
args = commandArgs(trailingOnly=TRUE)
species <- args[1]
refname <- args[2]
sampname <- args[3]

# load libraries
suppressMessages(library(dplyr))

# source file with all necessary functions
# make sure this file is in the same directory
source('./LR_functions.R')

# read in reference file and sample file
refgts <- read.table(refname, sep=',', header=TRUE)
sampgts <- read.table(sampfname, sep=',', header=TRUE)

# use same theta estimates as in the 2018 Sci Adv paper
if (Species == "forest"){
  Theta <- 0.059
} else{
  Theta <- 0.047
}

# generate the results
results <- run_calc(refgts, sampgts, theta=Theta)

# write the output
write.table(results, file=paste0("obsLRs.", Species, ".txt"), 
            sep = "\t", 
            quote = F, 
            row.names = F)
