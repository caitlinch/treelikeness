# R functions to run parametric bootstrap

# Given a folder that contains an IQ-tree output, extract the model used, alignment length and number of taxa
get.simulation.parameters <- function(folder_path){
  # Read in the summary information about the alignment
  alignment_info <- read.table(paste0(folder_path,"/alignment.nex-seq-summary.txt"),header = TRUE)
  # Extract the number of taxa 
  ntaxa <- nrow(alignment_info)
  # Extract the length of the alignment - take maximum length of all input sequences
  nsites <- max(alignment_info$Sequence_length)
  # read in the IQ-TREE file to get substitution model and parameters
  iqtree_file <- paste0(folder_path,"/alignment.nex.iqtree")
  iq_db <- readLines(iqtree_file)
  start_ind      <- which(iq_db == "SUBSTITUTION PROCESS")
  model_ind      <- start_ind+3
  rate_ind       <- start_ind+5
  AC_ind         <- start_ind+7
  AG_ind         <- start_ind+8
  AT_ind         <- start_ind+9
  CG_ind         <- start_ind+10
  CT_ind         <- start_ind+11
  GT_ind         <- start_ind+12
  a_freq         <- start_ind+16
  c_freq         <- start_ind+17
  g_freq         <- start_ind+18
  t_freq         <- start_ind+19
  rate_matrix    <- (start_ind+23):(start_ind+26)
  model_rate_het <- (start_ind+28):(start_ind+29)
  # cat matrix has variable number of rows, but starts 31 rows after start_ind
  cat_start <- start_ind+31
  # first find all places where there's an empty row
  empty_rows <- which(iq_db == "")
  empty_rows <- empty_rows[empty_rows > cat_start]
  # find closest number to start of matrix: this will be the empty row after the matrix finishes
  cat_end <- which.min(abs(empty_rows)-cat_start)
  cat_matrix <- cat_start:cat_end
  
  # parameters in order: number of taxa, number of sites, 
  pms <- list(ntaxa,nsites)
  return(pms)
}

# Given the relevant information, run one parametric bootstrap
do.1.bootstrap <- function(folder_path,parameters) {
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

  # 4. Calculate tree likeness of that alignment
    # need to read in nexus to IQ-tree?
    # to convert to nexus: write.nexus.data(dna_sim, file = "/Users/caitlin/Downloads/test.nexus")
}

