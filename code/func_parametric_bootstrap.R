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
    # scale tree to have total depth of 0.6 (CHECK THIS VALUE & FIND A BIOLOGICAL REASON FOR IT)
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
    rate1 <- as.numeric(strsplit(iq_file[[grep("A-C",iq_file)]],":")[[1]][2]) # A-C rate (same as code above, but combined 4 lines into 1 line)
    rate2 <- as.numeric(strsplit(iq_file[[grep("A-G",iq_file)]],":")[[1]][2]) # A-G rate
    rate3 <- as.numeric(strsplit(iq_file[[grep("A-T",iq_file)]],":")[[1]][2]) # A-T rate
    rate4 <- as.numeric(strsplit(iq_file[[grep("C-G",iq_file)]],":")[[1]][2]) # C-G rate
    rate5 <- as.numeric(strsplit(iq_file[[grep("C-T",iq_file)]],":")[[1]][2]) # C-T rate
    rate6 <- as.numeric(strsplit(iq_file[[grep("G-T",iq_file)]],":")[[1]][2]) # G-T rate
    
    # Extract the state frequencies
    sf1 <- as.numeric(strsplit(iq_file[[grep("pi\\(A\\)",iq_file)]],"=")[[1]][2]) # pi(A) - A freq. Remember to double backslash to escape before brackets
    sf2 <- as.numeric(strsplit(iq_file[[grep("pi\\(C\\)",iq_file)]],"=")[[1]][2]) # pi(C) - C freq
    sf3 <- as.numeric(strsplit(iq_file[[grep("pi\\(G\\)",iq_file)]],"=")[[1]][2]) # pi(G) - G freq
    sf4 <- as.numeric(strsplit(iq_file[[grep("pi\\(T\\)",iq_file)]],"=")[[1]][2]) # pi(T) - T freq
    
    # Extract model of rate heterogeneity
    mrh1      <- strsplit(iq_file[[grep("Model of rate heterogeneity:",iq_file)]],":")[[1]][2] # Extract model of rate heterogeneity 
    mrh2      <- as.numeric(strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][2]) # Line after the "model of rate heterogeneity" varies - extract it regardless of what it is 
    mrh2_name <- strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][1] # As the line varies, extract the name for the output dataframe
    mrh2_name <- gsub(" ","_",mrh2_name) # change the name to be easy to parse
    
    # make a list of the rows for the output dataframe
    names <- c("file_name","sequence_type","n_taxa","n_sites","substitution_model", "A-C_rate","A-G_rate","A-T_rate","C-G_rate","C-T_rate","G-T_rate","A_freq","C_freq",
               "G_freq","T_freq","model_of_rate_heterogeneity",mrh2_name)
    # Make a list of the output rows for the output dataframe
    op <- c(op1,"DNA",op2,op3,op4,rate1,rate2,rate3,rate4,rate5,rate6,sf1,sf2,sf3,sf4,mrh1,mrh2)
    # Create the output dataframe
    op_df <- data.frame(names,op)
    # Name the columns
    names(op_df) <- c("parameter","value")
    
    # Create the rate matrix Q
    Q_start <- grep("Rate matrix Q:",iq_file)+2
    Q_end   <- Q_start+3
    # Create the columns
    c1 <- c("A","C","G","T")
    c2 <- c()
    c3 <- c()
    c4 <- c()
    c5 <- c()
    # For each row in the iqtree file rate matrix
    for (i in Q_start:Q_end){
      # Split the row
      row <- strsplit(iq_file[[i]],"   ")
      # Add the resulting values to the relevant columns
      c2 <- c(c2,as.numeric(row[[1]][2])) # convert to numeric so can use the numbers more easily later
      c3 <- c(c3,as.numeric(row[[1]][3]))
      c4 <- c(c4,as.numeric(row[[1]][4]))
      c5 <- c(c5,as.numeric(row[[1]][5]))
    }
    # Create a dataframe of the rate matrix Q
    q_df <- data.frame(c1,c2,c3,c4,c5)
    #Rename the columns
    names(q_df) <- c("nucleotide","A","C","G","T")
    
    #Create the matrix for discrete gamma categories
    g_start <- grep(" Category",iq_file)+1 # get the index for the first line of the gamma categories matrix
    empty   <- which(iq_file=="") # get indexes of all empty lines
    empty   <- empty[empty>g_start] # get empty lines above gamma categories matrix
    g_end   <- empty[1]-1 # get end index for gamma categories matrix (one less than next empty line)
    g1 <- c() # initialise columns to store data in
    g2 <- c()
    g3 <- c()
    # Iterate through rows in gamma matrix
    for (i in g_start:g_end){
      row <- strsplit(iq_file[[i]],"        ") # split the rows on the long strong of 0's in the middle
      g1 <- c(g1,as.numeric(row[[1]][1])) # add the values to the columns
      g2 <- c(g2,as.numeric(row[[1]][2]))
      g3 <- c(g3,as.numeric(row[[1]][3]))
    }
    g_df <- data.frame(g1,g2,g3) # create a dataframe of the information
    names(g_df) <- c("category","relative_rate","proportion") # name the columns
    
    # Create a list of the three dataframes
    # This will be the output 
    params <- list(op_df,g_df,q_df)
    # Name the parameters so they're easy to access once you've outputted the data
    names(params) <- c("parameters","gamma_categories","Q_rate_matrix")
    
  } else if (input_ls[[1]][7]=="amino-acid"){ # alternatively if the data is amino acid sites
    # Extract model of rate heterogeneity
    mrh1      <- strsplit(iq_file[[grep("Model of rate heterogeneity:",iq_file)]],":")[[1]][2] # Extract model of rate heterogeneity 
    mrh2      <- as.numeric(strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][2]) # Line after the "model of rate heterogeneity" varies - extract it regardless of what it is 
    mrh2_name <- strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][1] # As the line varies, extract the name for the output dataframe
    mrh2_name <- gsub(" ","_",mrh2_name) # change the name to be easy to parse
    # Extract state frequencies
    sf1      <- strsplit(iq_file[[grep("State frequencies:",iq_file)]],":")[[1]][2]
    
    # make a list of the rows for the output dataframe
    names <- c("file_name","sequence_type","n_taxa","n_sites","substitution_model","model_of_rate_heterogeneity",mrh2_name,"state_frequencies")
    # Make a list of the output rows for the first output dataframe
    op <- c(op1,"amino-acid",op2,op3,op4,mrh1,mrh2,sf1)
    op_df <- data.frame(names,op)
    names(op_df) <- c("parameter","value")
    
    #Create the matrix for discrete gamma categories
    g_start <- grep(" Category",iq_file)+1 # get the index for the first line of the gamma categories matrix
    empty   <- which(iq_file=="") # get indexes of all empty lines
    empty   <- empty[empty>g_start] # get empty lines above gamma categories matrix
    g_end   <- empty[1]-1 # get end index for gamma categories matrix (one less than next empty line)
    g1 <- c() # initialise columns to store data in
    g2 <- c()
    g3 <- c()
    # Iterate through rows in gamma matrix
    for (i in g_start:g_end){
      row <- strsplit(iq_file[[i]],"        ") # split the rows on the long strong of 0's in the middle
      g1 <- c(g1,as.numeric(row[[1]][1])) # add the values to the columns
      g2 <- c(g2,as.numeric(row[[1]][2]))
      g3 <- c(g3,as.numeric(row[[1]][3]))
    }
    g_df <- data.frame(g1,g2,g3) # create a dataframe of the information
    names(g_df) <- c("category","relative_rate","proportion") # name the columns
    
    # Create a list of the dataframes - this will be the output
    params <- list(op_df,g_df)
    # Name the parameters si tget're easy to access once you've outputted the data
    names(params) <- c("parameters","gamma_categories")
  }
  
  # Now the information has been collected, create an output dataframe
  # make the output dataframe
  return(params)
}

