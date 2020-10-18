#!/usr/local/bin/Rscript

## load user-defined functions, packages
suppressMessages(library("argparse"))
suppressMessages(library("dplyr"))
source("./do_one.R")
source("./LR_functions.R")

## load command line arguments
parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "elephants",
                    help = "name of simulation")
## both of the below should be 1 at all times for this project
parser$add_argument("--nreps-total", type = "double", default = 1,
                    help = "number of replicates for each set of params")
parser$add_argument("--nreps-per-job", type = "double", default = 1,
                    help = "number of replicates per job")
args <- parser$parse_args()

## set up a grid of parameters to cycle over
# number of loci
ns <- seq(10,16)
# type of LR calculated
types <- c("FS", "HS", "PO", "DM")
# species
specs <- c("forest", "savannah")

# number of jobs for each parameter combo
njobs_per_combo <- args$nreps_total/args$nreps_per_job
# set up grid of parameters
param_grid <- expand.grid(mc_id = 1:njobs_per_combo, 
                          nloci = ns,
                          type = types,
                          species = specs)

## get job id from scheduler
job_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
## current dynamic args
current_dynamic_args <- param_grid[job_id, ]

## set seed
current_seed <- job_id + current_dynamic_args$n
set.seed(current_seed)

## each simulation setting generates 1 million LRs
niters <- 1000000
for (i in 1:args$nreps_per_job){
  do_one(nloci = current_dynamic_args$nloci,
         niters = niters,
         type = current_dynamic_args$type,
         species = current_dynamic_args$species)
}
