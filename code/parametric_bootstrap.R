# R functions to run parametric bootstrap

# Given the relevant information, run one parametric bootstrap
do.1.bootstrap <- function(iq_path,folder_path,parameters,test_statistic) {
  params <- get.simulation.parameters(folder_path)
  ntaxa <- 
  nsites <-
  model <- 
  
  # Start the parametric bootstrap process
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
  # save alignment as a nexus file
  op_file <- paste0(folder_path,"/alignment.nexus")
  write.nexus.data(dna_sim, file = op_file)
  
  # 4. Calculate tree likeness of that alignment
  # calculate the test statistic using numbers from the grant proposal  
  if (test_statistic ==1){
    ts <- pdm.ratio(iq_path,op_file)
  } else if (test_statistic == 2){
    ts <- normalised.pdm.ratio(iq_path,op_file)
  } else if (test_statistic == 3){
    print("haven't written the third test statistic yet")
  } else{
    break
  }
}

# Given a .iqtree file, a number of replicates and an output directory, this function will
# generate that number of alignments in the output directory with the parameters from 
# the .iqtree file (number of taxa, number of sites, rates
# and base frequencies)



# Given a .iqtree file, this function will extract the relevant parameters
get.simulation.parameters <- function(folder_path){
  # Read in the summary information about the alignment
  alignment_info <- read.table(paste0(folder_path,"/alignment.nex-seq-summary.txt"),header = TRUE)
  # Extract the number of taxa 
  ntaxa <- nrow(alignment_info)
  # Extract the length of the alignment - take maximum length of all input sequences
  nsites <- max(alignment_info$Sequence_length)
  # read in the IQ-TREE file to get substitution model and parameters
  iqtree_file <- paste0(folder_path,"/alignment.nex.iqtree")
  iq_df <- readLines(iqtree_file)
  # get the start and end of the substitution parameters section in the iqtree file
  start_ind      <- which(iq_df == "SUBSTITUTION PROCESS")+3
  end_ind        <- which(iq_df == "MAXIMUM LIKELIHOOD TREE")-2
  # extract the substitution model lines of the iq_df
  subs_df <- iq_df[start_ind:end_ind]
  # parameters in order: number of taxa, number of sites, substitution database
  pms <- list(ntaxa,nsites,subs_df)
  return(pms)
}

