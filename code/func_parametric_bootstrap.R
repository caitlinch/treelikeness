# R functions to run parametric bootstrap

library(phytools)
library(seqinr)
library(ape)
library(phangorn)

# Given a directory and a number of replicates, this function will
# generate that number of alignments in the directory with the parameters from 
# the .iqtree file (number of taxa, number of sites, rates
# and base frequencies)
# make sure the alignment folder ends with a slash!
phylo.bootstrap <- function(alignment_folder,n_reps,iq_path,splitstree_path) {
  # extract the name of the .iqtree file that contains the parameters for the simulation
  dotiqtree_path <- paste0(folder_path, list.files(folder_path)[grep("iqtree", list.files(folder_path))])
  params <- get.simulation.parameters(dotiqtree_path) #need to feed in .iqtree file
  
  # fill out the rep numbers (padded with 0s to get to 4 digits)
  rep_ids <- 1:n_reps
  rep_ids <- sprintf("%04d", rep_ids)

}

# Given the relevant information, run one parametric bootstrap (create the alignment and run the test statistics, output the p-values as a vector)
do.1.bootstrap <- function(rep_number,params,alignment_folder,iq_path,splitstree_path){
  # Create a new folder name to store this alignment and its outputs in
  bs_folder <- paste0(alignment_folder,rep_number,"/")
  
  if (dir.exists(bs_folder)==TRUE){
    redo <- TRUE
  } else if (dir.exists(bs_folder)==FALSE){
    # if the folder doesn't exist, create it
    dir.create(bs_folder)
    redo <- FALSE  # if the directory and the alignment don't exist, this is not a redo and the alignment needs to be tested
  } 
  
  # Only create alignments, run IQ-Tree and all the tests/test statistics and everything if this is NOT a redo (redo == FALSE)
  if (redo == FALSE){
    # Create the alignment 
    # Sample code for generating a parametric DNA sequence if you have a tree
    # s1 = simSeq(t1, l = 500, type="DNA", bf=c(.25,.25,.25,.25), Q=c(1,1,1,1,1,1), rate=1)
    # s2 = simSeq(t2, l = 500, type="DNA", bf=c(.25,.25,.25,.25), Q=c(1,1,1,1,1,1), rate=1)
    # aln = c(s1, s2)
    
    # Extract the parameters you need to enter into simSeq
    n_bp = 1300 # want a sequence that is 1300 base pairs long
    # Extract the vector form of the rate matrix 
    m <- params$Q_rate_matrix[,2:5] # extract the square block with the rates and not the header column
    Q_vec <- c(m[2,1],m[3,1],m[3,2],m[4,1],m[4,2],m[4,3]) # extract the rates 
    # Extract the base frequencies in the following order: A, C, G, T
    base_freqs <- c(as.numeric(params$parameters[[12,2]]), as.numeric(params$parameters[[13,2]]), as.numeric(params$parameters[[14,2]]), as.numeric(params$parameters[[15,2]]))
    seq_type <- "DNA" # generate DNA sequence
      
    # The alignment now definitely exists. Now you can run IQ-tree on the alignment
    call.IQTREE.quartet(program_paths[["IQTree"]],al_file,row[["n_taxa"]])
  }
  
  # Calculate the test statistics
}


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
    state_freq_line <- iq_file[[grep("State frequencies",iq_file)]]
    if (state_freq_line == "State frequencies: (equal frequencies)"){
      # If the state frequencies are all equal, assign them all to 0.25 (1/4)
      sf1 <- 0.25 # pi(A) - A freq.
      sf2 <- 0.25 # pi(C) - C freq.
      sf3 <- 0.25 # pi(G) - G freq.
      sf4 <- 0.25 # pi(T) - T freq.
    } else {
      # If the state frequencies are not all equal, extract what they are
      sf1 <- as.numeric(strsplit(iq_file[[grep("pi\\(A\\)",iq_file)]],"=")[[1]][2]) # pi(A) - A freq. Remember to double backslash to escape before brackets
      sf2 <- as.numeric(strsplit(iq_file[[grep("pi\\(C\\)",iq_file)]],"=")[[1]][2]) # pi(C) - C freq.
      sf3 <- as.numeric(strsplit(iq_file[[grep("pi\\(G\\)",iq_file)]],"=")[[1]][2]) # pi(G) - G freq.
      sf4 <- as.numeric(strsplit(iq_file[[grep("pi\\(T\\)",iq_file)]],"=")[[1]][2]) # pi(T) - T freq.
    }
    
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
    op_df <- data.frame(names,op, stringsAsFactors = FALSE)
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
      row <- strsplit(iq_file[[i]],"   ")[[1]]
      row <- row[row!=""] # remove any empty strings from the vector
      row <- row[row!=" "] # remove any 1 space strings from the vector
      # Add the resulting values to the relevant columns
      c2 <- c(c2,as.numeric(row[2])) # convert to numeric so can use the numbers more easily later
      c3 <- c(c3,as.numeric(row[3]))
      c4 <- c(c4,as.numeric(row[4]))
      c5 <- c(c5,as.numeric(row[5]))
    }
    # Create a dataframe of the rate matrix Q
    q_df <- data.frame(c1,c2,c3,c4,c5, stringsAsFactors = FALSE)
    #Rename the columns
    names(q_df) <- c("nucleotide","A","C","G","T")
    
    # Check if the model for rate heterogeneity is uniform
    mrh1_check <- gsub(" ","",mrh1)
    if (mrh1_check=="Uniform"){
      # If the model for rate heterogeneity is uniform, don't need to create a matrix for discrete gamma rate categories
      g_df <- "Uniform"
    } else {
      # If the model isn't uniform, need to create a matrix to collect and store the gamme category information
      #Create the matrix for discrete gamma categories
      g_start <- grep(" Category",iq_file)+1 # get the index for the first line of the gamma categories matrix
      empty   <- which(iq_file=="") # get indexes of all empty lines
      empty   <- empty[empty>g_start] # get empty lines above gamma categories matrix
      g_end   <- empty[1]-1 # get end index for gamma categories matrix (one less than next empty line)
      end_line <- iq_file[g_end]
      # if the end isn't an empty line, subtract one from the end count 
      # to exclude lines like "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category."
      if (end_line != ""){
        g_end = g_end - 1
      }
      # Start collecting info for the matrix
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
      g_df <- data.frame(g1,g2,g3, stringsAsFactors = FALSE) # create a dataframe of the information
      names(g_df) <- c("category","relative_rate","proportion") # name the columns
    }
    
    
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

