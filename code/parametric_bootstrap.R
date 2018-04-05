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
get.simulation.parameters <- function(dotiqtree_file){
  # read in the IQ-TREE file to get substitution model and parameters
  iq_file <- readLines(dotiqtree_file)
  # extract the file name
  ind      <- grep("Input file name:",iq_file)
  op1      <- substr(iq_file[[ind]],18,nchar(iq_file[[ind]]))
  # extract the number of taxa and extract the length of the alignment
  ind         <- grep("Input data:",iq_file)
  input_str   <- iq_file[[ind]] # get the line that contains this info
  input_ls    <- strsplit(input_str," ")
  op2         <- input_ls[[1]][3] # extract number of sequences (number of taxa)
  op3         <- input_ls[[1]][6] # extract number of sites 
  # Extract the model of substitution (same for amino-acid and nucleotide files)
  ind         <- grep("Model of substitution: ",iq_file)
  sub_str     <- iq_file[[ind]] # get the line that contains this info
  sub_ls      <- strsplit(sub_str," ")
  op4         <- sub_ls[[1]][4]
  
  # check the type of sites - amino acid or DNA
  # if the input is DNA (nucleotide sites), gather that information
  if (input_ls[[1]][7]=="nucleotide"){
    # Extract the rate parameters
    # A-C: 1.1021
    # A-G: 3.5984
    # A-T: 0.6766
    # C-G: 1.4170
    # C-T: 3.5984
    # G-T: 1.0000
  
    
    # make a list of the rows for the output dataframe
    names <- c("file_name","n_taxa","n_sites","substitution_model", "A-C_rate","A-G_rate","A-T_rate","C-G_rate","C-T_rate","G-T_rate")
    # Make a list of the output rows for the output dataframe
    op <- c(op1,op2,op3,op4)
  } else if (input_ls[[1]][7]=="amino-acid"){ # alternatively if the data is amino acid sites
    # do nothing for now whoops haha (add this in after you've done the DNA part)
    # make a list of the rows for the output dataframe
    names <- c("file_name","n_taxa","n_sites","substitution_model")
    # Make a list of the output rows for the output dataframe
    op <- c(op1,op2,op3,op4)
  }
  
  # Now the information has been collected, create an output dataframe
  # make the output dataframe
  op_df <- data.frame(names,op)
  names(op_df) <- c("parameter","value")
  return(op_df)
}

