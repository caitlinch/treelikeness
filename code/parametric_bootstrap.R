parametric.bootstrap1 <- function(ntaxa,nsites = 1000){
  # 1. Simulate a tree
    # simulate a birth-death tree on a fixed number of extant taxa
        # n = number of taxa, numbsim = # of trees to simulate, lambda = speciation rate
        # mu = extinction rate, 
  tree_sim <- sim.bd.taxa(ntaxa, numbsim = 1, lambda = 0.5, mu = 0, frac = 1, complete = FALSE, stochsampling = TRUE)[[1]]
  
  # 2. Simulate rate variation
    # set parameters for creating rate variation in tree
  mol_rate <- 0.1
  mrate_sd <- 0.1
    # make the rate variation in the tree
  rate_var <- tree_sim
  rate_var$edge.length <- tree_sim$edge.length * rlnorm(length(tree_sim$edge.length), meanlog = log(mol_rate), sdlog = mrate_sd)
    # scale tree to have total depth of 0.6 million years (CHECK THIS VALUE & FIND A BIOLOGICAL REASON FOR IT)
  rate_var <- rescale(rate_var,"depth",0.6)

  # 3. Simulate an alignment 
  if (nsites > 0){
    # simulate the sequence alignment based on the tree and number of sites
    dna_sim <- as.DNAbin(simSeq(rate_var, l = nsites))
  } else {
    # if there are no sites, create a dummy matrix 
    dna_sim <- matrix('n', ncol = 1, nrow = n_taxa)
    rownames(dna_sim) <- tree_sim$tip.label
    dna_sim <- as.DNAbin(dna_sim)
  }

  # 4. Calculate tree likeness of that alignment
    # need to read in nexus to IQ-tree?
    # to convert to nexus: write.nexus.data(dna_sim, file = "/Users/caitlin/Downloads/test.nexus")
}

library(TreeSim)
library(ape)
library(phytools)
library(phangorn)
ntaxa <- 10
