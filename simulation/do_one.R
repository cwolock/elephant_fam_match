# run 1 simulation setting
do_one <- function(nloci, niters, type, species){
  fname <- paste0("elephant_ref_", 
                  species, 
                  ".filtered.longformat.", 
                  nloci, 
                  ".csv") 
  gts <- read.table(fname, sep=',', header=TRUE)
  # use theta estimates from the 2018 Sciences Advances paper
  if (species == "forest"){
    theta <- 0.059
  } else{
    theta <- 0.047
  }
  results <- run_sim(gts, theta=theta, niters=niters, 
                     nloci=nloci, type=type)
  type = as.character(type)
  write.csv(results, file=paste0("LRs.", species, ".", 
                                  type, ".", niters, ".", 
                                  nloci, ".csv"))
}
