# set kinship kappas
UN_k <- c(1,0,0)
DM_k <- c(0,0,1)
PO_k <- c(0,1,0)
FS_k <- c(.25, .5, .25)
HS_k <- c(0.5, 0.5, 0)

### calculates allele frequencies for a marker
calc_afs <- function(genotypes, allele_names){
  # drop -999 which corresponds to missing
  genotypes <- genotypes[!is.na(genotypes)]
  afs <- rep(NA, length(allele_names))
  for (i in 1:length(allele_names)){
    freq <- sum(genotypes == allele_names[i]) / length(genotypes)
    afs[i] <- freq
  }
  return(afs)
}

### takes 2x2 matrix of parental alleles and samples one from each parent
### to create two "children"
inherit <- function(parents){
  kid1 <- sample(parents[1,1:2], 1, prob=c(.5, .5))
  kid2 <- sample(parents[2,1:2], 1, prob=c(.5, .5))
  return(c(kid1,kid2))
}

### calculates LR, takes 2x2 matrix of alleles
### we'll call the first row of the matrix G1, and the second row G2
calc_kinship <- function(unorder_gt_pair, kap, alleles, theta, afs){
  gt_pair <- matrix(NA, nrow=2, ncol=3)
  if (unorder_gt_pair[2,1] == unorder_gt_pair[2,2]){
    gt_pair[1,] <- unorder_gt_pair[2,]
    gt_pair[2,] <- unorder_gt_pair[1,]
  } else{
    gt_pair <- unorder_gt_pair
  }
  
  allele_name <- gt_pair[1,3]
  
  gt_pair <- gt_pair[,-3]
  
  gt_pair <- apply(gt_pair, 2, as.numeric)
  
  # get AFs from AF data frame
  p.a <- afs[which(alleles == gt_pair[1,1]),allele_name]
  p.b <- afs[which(alleles == gt_pair[1,2]),allele_name]
  p.c <- afs[which(alleles == gt_pair[2,1]),allele_name]
  p.d <- afs[which(alleles == gt_pair[2,2]),allele_name]
  # g1 homozygous
  if (gt_pair[1,1] == gt_pair[1,2]){
    # g1 = g2, both homozygous (CASE 1)
    if (sum(gt_pair[2,] == gt_pair[1,]) == 2){
      p_RE <- kap[1]*(p.a*(theta + (1-theta)*p.a)*(2*theta + 
              (1-theta)*p.a)*(3*theta + (1-theta)*p.a) / 
                          ((1+theta)*(1+2*theta))) + 
              kap[2]*(p.a*(theta + (1-theta)*p.a)*(2*theta + 
              (1-theta)*p.a) / (1+theta)) + 
              kap[3]*(p.a*(theta + (1-theta)*p.a))
      p_UN <- UN_k[1]*(p.a*(theta + (1-theta)*p.a)*(2*theta + 
              (1-theta)*p.a)*(3*theta + (1-theta)*p.a) / 
                         ((1+theta)*(1+2*theta))) + 
              UN_k[2]*(p.a*(theta + (1-theta)*p.a)*(2*theta + 
              (1-theta)*p.a) / (1+theta)) + 
              UN_k[3]*(p.a*(theta + (1-theta)*p.a))
    } else if(gt_pair[2,1] == gt_pair[1,1] | 
              gt_pair[2,2] == gt_pair[1,1]){
    # one matching allele (CASE 3)
      # find which allele from gt2 matched the pair from gt1
      if (gt_pair[2,1] == gt_pair[1,1]){
        p.diff <- p.d
      } else{
        p.diff <- p.c
      }
      p_RE <- 2*kap[1]*(p.a*(theta + (1-theta)*p.a)*(2*theta + 
              (1-theta)*p.a)*((1-theta)*p.diff) / 
                         ((1+theta)*(1+2*theta))) + 
              kap[2]*(p.a*(theta + (1-theta)*p.a)*((1-theta)*p.diff) / 
              (1+theta))
        
      p_UN <- 2*UN_k[1]*(p.a*(theta + (1-theta)*p.a)*(2*theta + 
              (1-theta)*p.a)*((1-theta)*p.diff) / 
                           ((1+theta)*(1+2*theta))) + 
              UN_k[2]*(p.a*(theta + (1-theta)*p.a)*((1-theta)*p.diff) / 
              (1+theta))
    } else{
    # no matching alleles
      # g2 homozygous (CASE 2)
      if (gt_pair[2,1] == gt_pair[2,2]){
        p_RE <- kap[1]*(p.a*(theta + (1-theta)*p.a)*
                (1-theta)*p.c*(theta + (1-theta)*p.c) / 
                            ((1+theta)*(1+2*theta)))
        p_UN <- UN_k[1]*(p.a*(theta + (1-theta)*p.a)*
                (1-theta)*p.c*(theta + (1-theta)*p.c) / 
                            ((1+theta)*(1+2*theta)))
      } else{
      # g2 heterozygous (CASE 4)
        p_RE <- 2*kap[1]*(p.a*(theta + (1-theta)*p.a)*(1-theta)*
                p.c*((1-theta)*p.d) / ((1+theta)*(1+2*theta)))
        p_UN <- 2*UN_k[1]*(p.a*(theta + (1-theta)*p.a)*(1-theta)*
                p.c*((1-theta)*p.d) / ((1+theta)*(1+2*theta)))
      }
    }
  } else{
  # g1 heterzygous
    # g1 = g2 (CASE 5)
    if (sum(gt_pair[2,] == gt_pair[1,]) == 2 | 
        sum(gt_pair[2,] == c(gt_pair[1,2], gt_pair[1,1])) == 2){
      p_RE <- 4*kap[1]*(p.a*(1-theta)*p.b*(theta + 
              (1-theta)*p.a)*(theta + (1-theta)*p.b) / 
              ((1+theta)*(1+2*theta))) + 
              kap[2]*(p.a*(1-theta)*p.b*((theta + (1-theta)*p.a) + 
              (theta + (1-theta)*p.b)) / (1+theta)) + 
              2*kap[3]*p.a*(1-theta)*p.b
      p_UN <- 4*UN_k[1]*(p.a*(1-theta)*p.b*(theta + 
              (1-theta)*p.a)*(theta + (1-theta)*p.b) / 
              ((1+theta)*(1+2*theta))) + 
              UN_k[2]*(p.a*(1-theta)*p.b*((theta + (1-theta)*p.a) + 
              (theta + (1-theta)*p.b)) / (1+theta)) + 
              2*UN_k[3]*p.a*(1-theta)*p.b
    } else if(gt_pair[2,1] == gt_pair[1,1] | 
              gt_pair[2,1] == gt_pair[1,2] | 
            gt_pair[2,2] == gt_pair[1,1] | 
            gt_pair[2,2] == gt_pair[1,2]){
    # g1 and g2 heterozygous with one shared allele (CASE 6)
      # find which alleles matched
      if (gt_pair[2,1] == gt_pair[1,1]){
        p.same <- p.a
        p.diffg1 <- p.b
        p.diffg2 <- p.d
      } else if(gt_pair[2,1] == gt_pair[1,2]){
        p.same <- p.b
        p.diffg1 <- p.a
        p.diffg2 <- p.d
      } else if(gt_pair[2,2] == gt_pair[1,1]){
        p.same <- p.a
        p.diffg1 <- p.b
        p.diffg2 <- p.c
      } else{
        p.same <- p.b
        p.diffg1 <- p.a
        p.diffg2 <- p.c
      }
      p_RE <- 4*kap[1]*(p.same*(1-theta)*p.diffg1*(theta + 
              (1-theta)*p.same)*((1-theta)*p.diffg2) / 
                          ((1+theta)*(1+2*theta))) + 
              kap[2]*(p.a*(1-theta)*p.diffg1*(1-theta)*p.diffg2 / 
              (1+theta))
      p_UN <- 4*UN_k[1]*(p.same*(1-theta)*p.diffg1*(theta + 
              (1-theta)*p.same)*((1-theta)*p.diffg2) / 
                          ((1+theta)*(1+2*theta))) + 
              UN_k[2]*(p.a*(1-theta)*p.diffg1*(1-theta)*p.diffg2 / 
              (1+theta))
    } else{
    # all four alleles different (CASE 7)
      p_RE <- 4*kap[1]*(p.a*(1-theta)*p.b*(1-theta)*p.c*(1-theta)*p.d / 
              ((1+theta)*(1+2*theta)))
      p_UN <- 4*UN_k[1]*(p.a*(1-theta)*p.b*(1-theta)*p.c*(1-theta)*p.d / 
              ((1+theta)*(1+2*theta)))
    }
  }
  LR <- p_RE / p_UN
  return(LR)
}

### function to run the simulations
run_sim <- function(pop_df, theta, niters, nloci, type){
  type = as.character(type)  
  pop_df[pop_df == -999] <- NA
  # get allele names (dropping -999 which corresponds to missing)
  alleles <- sort(unique(as.vector(t(data.matrix(pop_df)))))
  alleles <- alleles[!is.na(alleles)]
  # get a matrix of allele frequencies
  af_matrix <- apply(X = pop_df, MARGIN = 2, FUN = calc_afs, alleles)
  afs <- data.frame(af_matrix)
  
  # LRs for truly related pairs
  true_LRs <- rep(NA, niters)
  # LRs for randomly matched people
  false_LRs <- rep(NA, niters)
  # pick alleles of desired size, find 3 parents whohave those alleles
  for (i in 1:niters){
    # sample some loci
    index_to_include <- sample(1:ncol(pop_df), nloci, replace=FALSE)
    # filter data frame to only include those loci
    pop_df_copy <- pop_df[,index_to_include,drop=FALSE]
    # get rid of people with missing genotypes at the included loci
    pop_df_copy <- pop_df_copy[complete.cases(pop_df_copy),,drop=FALSE]
    # create empty family genotypes matrix (3 parents, 3 children)
    fam <- array(, c(length(index_to_include), 6, 3))
    fam[,,3] <- colnames(pop_df_copy)
    # get number of individuals
    npeople <- nrow(pop_df_copy)/2
    # select parent genotypes
    pids <- sample(seq(1, npeople), 3, replace=FALSE)
    prows <- c(2*pids[1]-1, 2*pids[1], 
               2*pids[2]-1, 2*pids[2], 
               2*pids[3]-1, 2*pids[3]) 
    pgts <- pop_df_copy[prows,,drop=FALSE]
    # fill in three parents genotypes
    fam[,1,1:2] <- t(pgts[1:2,,drop=FALSE]) 
    fam[,2,1:2] <- t(pgts[3:4,,drop=FALSE])
    fam[,3,1:2] <- t(pgts[5:6,,drop=FALSE])
    # select children genotypes from parents
    # 4/5 are FS (parents are 1/2), 5/6 are HS (6's parents are 2/3)
    fam[,4,1:2] <- t(apply(X = fam[,1:2,,drop=FALSE], 
                     MARGIN = 1, 
                     FUN = inherit))
    fam[,5,1:2] <- t(apply(X = fam[,1:2,,drop=FALSE], 
                     MARGIN = 1, 
                     FUN = inherit))
    fam[,6,1:2] <- t(apply(X = fam[,2:3,,drop=FALSE], 
                     MARGIN = 1, 
                     FUN = inherit))
    
    # randomly sample parent and offspring
    parent_ind <- sample(c(1,2), 1)
    offspring_ind <- sample(c(4,5), 1)
    # calculate the parent-offspring and sib LRs
    locus_FS_LRs <- apply(X = fam[,4:5,,drop=FALSE], 
                          MARGIN = 1,
                          FUN = calc_kinship, 
                          FS_k, alleles, theta, afs)
    locus_HS_LRs <- apply(X = fam[,5:6,,drop=FALSE], 
                          MARGIN = 1, 
                          FUN = calc_kinship, 
                          HS_k, alleles, theta, afs)
    locus_PO_LRs <- apply(X = fam[,c(parent_ind, offspring_ind),,
                                   drop=FALSE], 
                          MARGIN = 1, 
                          FUN = calc_kinship, 
                          PO_k, alleles, theta, afs)
    # replace parent 2 by identical parent 1 for direct match
    fam[,2,] <- fam[,1,]
    locus_DM_LRs <- apply(X = fam[,1:2,,drop=FALSE], 
                    MARGIN = 1, 
                    FUN = calc_kinship, 
                    DM_k, alleles, theta, afs)
    
    FS_LR <- prod(locus_FS_LRs)
    HS_LR <- prod(locus_HS_LRs)
    PO_LR <- prod(locus_PO_LRs)
    DM_LR <- prod(locus_DM_LRs)
    
    if (type == "PO"){
      true_LRs[i] <- PO_LR
    } else if (type == "FS"){
      true_LRs[i] <- FS_LR
    } else if (type == "HS"){
      true_LRs[i] <- HS_LR
    } else{
      true_LRs[i] <- DM_LR
    }
    
    # now randomly select an unrelated pair
    pair <- array(, c(length(index_to_include), 2, 3))
    pair[,,3] <- colnames(pop_df_copy)
    ids <- sample(seq(1, npeople), 
                  2, 
                  replace=FALSE) # select two person ids
    rows <- c(2*ids[1]-1, 2*ids[1], 2*ids[2]-1, 2*ids[2])
    pgts <- pop_df_copy[rows,,drop=FALSE]
    pair[,1,1:2] <- t(pgts[1:2,,drop=FALSE]) 
    pair[,2,1:2] <- t(pgts[3:4,,drop=FALSE])
        
    # calculate the LRs for UN individuals
    locus_FS_LRs <- apply(X = pair, MARGIN = 1, 
                          FUN = calc_kinship, 
                          FS_k, alleles, theta, afs) 
    locus_HS_LRs <- apply(X = pair, MARGIN = 1, 
                          FUN = calc_kinship, 
                          HS_k, alleles, theta, afs) 
    locus_PO_LRs <- apply(X = pair, MARGIN = 1, 
                          FUN = calc_kinship, 
                          PO_k, alleles, theta, afs)
    locus_DM_LRs <- apply(X = pair, MARGIN = 1, 
                          FUN = calc_kinship, 
                          DM_k, alleles, theta, afs)
    FS_LR <- prod(locus_FS_LRs)
    HS_LR <- prod(locus_HS_LRs)
    PO_LR <- prod(locus_PO_LRs)
    DM_LR <- prod(locus_DM_LRs)
    if (type == "PO"){
      false_LRs[i] <- PO_LR
    } else if (type == "FS"){
      false_LRs[i] <- FS_LR
    } else if (type == "HS"){
      false_LRs[i] <- HS_LR
    } else{
      false_LRs[i] <- DM_LR
    }
  }
  # organize output
  dat <- data.frame(values = c(true_LRs, false_LRs),
                    type = c(rep(type, length(true_LRs)), 
                    rep('UN', length(false_LRs))))
  
  return(dat)
}
