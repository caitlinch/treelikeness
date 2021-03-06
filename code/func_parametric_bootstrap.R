# R functions to run parametric bootstrap

library(phytools)
library(seqinr)
library(ape)
library(phangorn)
library(stringr)

# Given a directory and a number of replicates, this function will
# generate that number of alignments in the directory with the parameters from 
# the .iqtree file (number of taxa, number of sites, rates
# and base frequencies)
# make sure the alignment folder ends with a slash!
phylo.parametric.bootstrap <- function(alignment_folder,n_reps,iq_path,splitstree_path, phipack_path, threeseq_path,exec_paths,tree_folder) {
  # extract the name of the .iqtree file that contains the parameters for the simulation
  dotiqtree_path <- paste0(alignment_folder, list.files(alignment_folder)[grep("iqtree", list.files(alignment_folder))])
  params <- get.simulation.parameters(dotiqtree_path) #need to feed in .iqtree file
  
  # Open the ML tree from IQ-Tree
  # extract the name of the file that contains the ML tree calculated by IQ-Tree for the original alignment (the alignment with recombination)
  all_tree_paths <- c(paste0(alignment_folder,"tree1.treefile"),paste0(alignment_folder,"tree2.treefile"),paste0(alignment_folder,"alignment.nexus.treefile"))
  ML_tree_path <- all_tree_paths[grep("alignment", all_tree_paths)] # get only the treefile for the alignment
  print(ML_tree_path)
  ML_tree <- read.tree(ML_tree_path)
  
  # fill out the rep numbers (padded with 0s to get to 4 digits)
  rep_ids <- 1:n_reps
  rep_ids <- sprintf("%04d", rep_ids)
  # run each rep
  lapply(rep_ids, do.1.bootstrap, params, ML_tree, alignment_folder, iq_path, splitstree_path, phipack_path, threeseq_path)
  
  # collate the bootstrap data and calculate the p-values.
  phylo.collate.bootstrap(alignment_folder,exec_paths, tree_folder)
}





# Given the relevant information, run one parametric bootstrap (create the alignment and run the test statistics, output the p-values as a vector)
do.1.bootstrap <- function(rep_number,params,tree,alignment_folder,iq_path,splitstree_path, phipack_path, threeseq_path){
  # Create a new folder name to store this alignment and its outputs in
  bs_folder <- paste0(alignment_folder,"bootstrap_",rep_number,"/")
  
  if (dir.exists(bs_folder)==TRUE){
    if (file.exists(paste0(bs_folder,"bootstrap_testStatistics.csv")) == TRUE) {
      redo <- TRUE # if the test statistics file exists, this is a redo and doesn't need doing
    } else {
      # The folder exists but doesn't have the bootstrap replicate test statistics - need to run those so make redo = FALSE (will not be redoing)
      redo = FALSE 
    }
  } else if (dir.exists(bs_folder)==FALSE){
    # if the folder doesn't exist, create it
    dir.create(bs_folder)
    redo <- FALSE  # if the directory and the alignment don't exist, this is not a redo and the alignment needs to be tested
  } 
  
  # Check if the alignment exists
  alignment.exists <- file.exists(paste0(bs_folder,"alignment.nexus"))
  
  # Change the working directory to the bootstrap folder
  setwd(bs_folder)
  # name the alignment 
  bs_al <- paste0(bs_folder, "alignment.nexus")
  
  if (alignment.exists == FALSE) {
    # If the alignment doesn't exist, create it

    # Sample code for generating a parametric DNA sequence if you have a tree
    # s1 = simSeq(t1, l = 500, type="DNA", bf=c(.25,.25,.25,.25), Q=c(1,1,1,1,1,1), rate=1)
    # s2 = simSeq(t2, l = 500, type="DNA", bf=c(.25,.25,.25,.25), Q=c(1,1,1,1,1,1), rate=1)
    # aln = c(s1, s2) # concatenate the alignments
    
    # Extract the parameters you need to enter into simSeq
    n_bp = as.numeric(params$parameters[4,2]) # want a sequence that is 1300 base pairs long
    # Extract the vector form of the rate matrix 
    m <- params$Q_rate_matrix[,2:5] # extract the square block with the rates and not the header column
    Q_vec <- c(m[2,1],m[3,1],m[3,2],m[4,1],m[4,2],m[4,3]) # extract the rates 
    # Extract the base frequencies in the following order: A, C, G, T
    base_freqs <- c(as.numeric(params$parameters[[12,2]]), as.numeric(params$parameters[[13,2]]), as.numeric(params$parameters[[14,2]]), as.numeric(params$parameters[[15,2]]))
    seq_type <- "DNA" # generate DNA sequence
    
    # Generate the DNA sequence
    # Don't need to specify rate variation because using JC model with no rate options
    dna_sim <- simSeq(tree, l = n_bp, type = seq_type, bf = base_freqs, Q = Q_vec)
    
    # Save the DNA alignment
    write.phyDat(dna_sim,file = bs_al, format = "nexus",interleaved = TRUE, datablock = FALSE) # write the output as a nexus file
    # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
    nexus <- readLines(bs_al) # open the new nexus file
    ind <- grep("BEGIN CHARACTERS",nexus)+2 # find which line
    nexus[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
    writeLines(nexus,bs_al) # output the edited nexus file
  }
  
  # Only create alignments, run IQ-Tree and all the tests/test statistics and everything if this is NOT a redo (redo == FALSE)
  if (redo == FALSE){
    # The alignment now definitely exists. Now you can run IQ-tree on the alignment
    n_taxa <- as.numeric(params$parameters[3,2]) # extract the number of taxa from the parameters 
    call.IQTREE.quartet.bootstrap(iq_path, bs_al, n_taxa)
    
    ## Calculate the test statistics
    # check is file is nexus or fasta
    filetype = tail(strsplit(bs_al,"\\.")[[1]],n=1) # extract file format
    if (filetype == "nexus"){
      # Phipack only reads in Phylip or fasta format - need to convert if the alignment is a nexus file
      data = read.nexus.data(bs_al) # read in nexus format alignment
      fasta.name <- paste0(bs_al,".fasta") # make a name for the fasta alignment by adding .fasta (super original ;) )
      write.fasta(sequences = data,names = names(data), file.out = fasta.name) # output alignment as a fasta format
      al_call_name <- fasta.name
    } else if (filetype == "fasta"){
      al_call_name <- bs_al
    }
    
    # run PHIPACK and 3seq (only if you haven't already - test by checking if the log files exist)
    phi_path <- phipack_path # get path to phipack executable
    seq_path <- threeseq_path # get path to 3seq executable
    
    if (file.exists(paste0(bs_folder,"3s.log")) == FALSE){
      seq_command <- paste0(seq_path," -f ", al_call_name)
      system(seq_command) #call 3SEQ
    }
    if (file.exists(paste0(bs_folder,"Phi.log")) == FALSE){
      phi_command <- paste0(phi_path," -f ",al_call_name, " -v") # assemble system command
      system(phi_command) #call phipack
    }
    
    # Extract significance from Phi Pack output
    phi_file <- paste0(bs_folder,"Phi.log")
    phi_file <- readLines(phi_file)
    ind      <- grep("p-Value",phi_file)
    phi_sig <- as.numeric(strsplit(phi_file[ind+3],":")[[1]][2])
    ind      <- grep("PHI Values",phi_file)
    phi_mean <- as.numeric(strsplit(phi_file[ind+4],"      ","")[[1]][2])
    phi_var <- as.numeric(strsplit(phi_file[ind+5],"      ","")[[1]][2])
    phi_obs <- as.numeric(strsplit(phi_file[ind+6],"      ","")[[1]][2])
    
    # Extract results output from 3Seq output
    seq_file <- paste0(bs_folder,"3s.log")
    seq_log <- readLines(seq_file) # open file
    ind      <- grep("Number of recombinant triplets",seq_log) # find the number of recombinant triplets line index
    num_trips <- seq_log[ind]
    num_trips <- strsplit(num_trips,":")[[1]][2] # extract the number of recombinant triplets
    num_trips <- trimws(num_trips) # trim the whitespace from the number of triplets
    ind      <- grep("Number of distinct recombinant sequences",seq_log) # find the number of distinct recombinant sequences line index
    num_dis <- seq_log[ind]
    num_dis <- strsplit(num_dis,":")[[1]][2] # extract the number of distinct recombinant sequences
    num_dis <- trimws(num_dis) # trim the whitespace from the number of distinct recombinant sequences
    # null hypothesis is of clonal evolution - need significant p-value to accept the alternative hypothesis
    ind      <- grep("Rejection of the null hypothesis of clonal evolution",seq_log) # find the p value line index
    seq_sig <- seq_log[ind]
    seq_sig <- strsplit(seq_sig,"=")[[1]][2] # extract the p value
    seq_sig <- trimws(seq_sig) # trim the whitespace from the number of distinct recombinant sequences
    
    # Extract quartet mapping from IQ-Tree results
    iq_log_path <- paste0(bs_al, ".iqtree")
    iq_log <- readLines(iq_log_path)
    ind <- grep("Number of fully resolved  quartets",iq_log)
    resolved_q <- as.numeric(strsplit(strsplit(iq_log[ind],":")[[1]][2],"\\(")[[1]][1])
    ind <- grep("Number of partly resolved quartets",iq_log)
    partly_resolved_q <- as.numeric(strsplit(strsplit(iq_log[ind],":")[[1]][2],"\\(")[[1]][1])
    ind <- grep("Number of unresolved",iq_log)
    unresolved_q <- as.numeric(strsplit(strsplit(iq_log[ind],":")[[1]][2],"\\(")[[1]][1])
    total_q <- (resolved_q+partly_resolved_q+unresolved_q)
    prop_resolved <- resolved_q/total_q
    
    # Calculate the median and mean delta score
    print(bs_al)
    pdmm <- as.matrix(mldist.pdm(bs_al)) # open pairwise distance matrix as a matrix
    deltaplot_results <- delta.plot(pdmm, k = 51, plot = FALSE) # calculate the delta.plot
    counts <- deltaplot_results$counts
    intervals <- seq(0,1,(1/(length(counts)-1)))
    deltaplot_df <- data.frame(intervals,counts)
    names(deltaplot_df) <- c("intervals","counts")
    deltaplot_df_name <- paste0(bs_folder,"deltaplot_histogram.csv")
    write.csv(deltaplot_df,file = deltaplot_df_name)
    # Want to calculate the mean and median delta q value - unfortunately the delta.plot function doesn't output raw data, so make a pseudo data set using the histogram values
    mean_dq <- mean(rep(deltaplot_df$intervals,deltaplot_df$counts)) # turn the interval data into a long list of "raw" values and calculate the mean
    median_dq <- median(rep(deltaplot_df$intervals,deltaplot_df$counts)) # turn the interval data into a long list of "raw" values and calculate the median
    mode_dq <- deltaplot_df[order(deltaplot_df$count, decreasing = TRUE),][1,1] # sort the dataframe by count values and extract the mode
    
    # Run my test statistics
    # Want the path to be the path to the alignment you're testing - here that's the bootstrap alignment (bs_al)
    # run pdm ratio (TS1) (modified splittable percentage)
    nn_trimmed <- tree.proportion(iqpath = iq_path, splitstree_path = splitstree_path, path = bs_al, network_algorithm = "neighbournet", trimmed = TRUE, seq_type = "dna")
    
    # Extract params csv from alignment folder
    all_files <- list.files(alignment_folder) # get a list of all the files
    ind <- grep("params",all_files) # find which of those files is the parameters file
    params_csv <- read.csv(paste0(alignment_folder,all_files[ind])) # open the parameter file
    params_csv <- params_csv[,2:ncol(params_csv)]
    
    # add test statistic results to params_csv
    params_csv$PHI_mean <- phi_mean
    params_csv$PHI_variance <- phi_var
    params_csv$PHI_observed <- phi_obs
    params_csv$X3SEQ_num_recombinant_triplets <- num_trips
    params_csv$X3SEQ_num_distinct_recombinant_sequences <- num_dis
    params_csv$LM_num_quartets <- total_q
    params_csv$LM_num_resolved_quartets <- resolved_q
    params_csv$LM_prop_resolved_quartets <- prop_resolved
    params_csv$LM_num_partially_resolved_quartets <- partly_resolved_q
    params_csv$LM_num_unresolved_quartets <- unresolved_q
    params_csv$tree_proportion <- nn_trimmed
    params_csv$mean_delta_q <- mean_dq
    params_csv$median_delta_q <- median_dq
    params_csv$mode_delta_q <- mode_dq
    params_csv$bootstrap_id <- paste0("bootstrap_",rep_number)
    
    # Output results as a csv in the bootstrap folder
    bs_results_csv <- paste0(bs_folder,"bootstrap_testStatistics.csv")
    write.csv(params_csv,file = bs_results_csv)
  }
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
  # Extract information about the sequence alignment
  ind <- grep("Number of constant sites:",iq_file)
  num_lines <- iq_file[c(ind:(ind+3))] # take the four lines listing the number of different kinds of sites
  num_split <- unlist(strsplit(num_lines,":")) # Split the lines at the colon
  num_names <- num_split[c(TRUE,FALSE)] # before the colon = name
  num_vals <- num_split[c(FALSE,TRUE)] # after the colon = value
  num_vals_regx <- regmatches(num_vals, gregexpr("[[:digit:]]+", num_vals)) # extract all the numbers after the colon
  # four strings = four lists of numbers: take first value from each list to get number of sites
  num_vals <- c(num_vals_regx[[1]][1],num_vals_regx[[2]][1],num_vals_regx[[3]][1],num_vals_regx[[4]][1]) 
  
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
    names <- c("file_name","sequence_type","n_taxa","n_sites",num_names,"substitution_model", "A-C_rate","A-G_rate","A-T_rate","C-G_rate","C-T_rate","G-T_rate","A_freq","C_freq",
               "G_freq","T_freq","model_of_rate_heterogeneity","model_of_rate_heterogeneity_line2_name","model_of_rate_heterogeneity_line2_value")
    # Make a list of the output rows for the output dataframe
    op <- c(op1,"DNA",op2,op3,num_vals,op4,rate1,rate2,rate3,rate4,rate5,rate6,sf1,sf2,sf3,sf4,mrh1,mrh2_name,mrh2)
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
      row <- strsplit(iq_file[[i]]," ")[[1]]
      row <- row[str_detect(row,"([0-9])")] # take only the numeric elements of the vector
      # Add the resulting values to the relevant columns
      c2 <- c(c2,as.numeric(row[1])) # convert to numeric so can use the numbers more easily later
      c3 <- c(c3,as.numeric(row[2]))
      c4 <- c(c4,as.numeric(row[3]))
      c5 <- c(c5,as.numeric(row[4]))
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
      # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
      check_line <- length(strsplit(strsplit(end_line, "        " )[[1]][1]," ")[[1]])
      if (check_line > 3){
        # If the check_line is longer than 3 characters, it won't be a group for the gamma categories but an instruction
        # Instructions can be excluded from the gamma matrix (but categories can't)
        g_end = g_end - 1
      }
      # Start collecting info for the matrix
      g1 <- c() # initialise columns to store data in
      g2 <- c()
      g3 <- c()
      # Iterate through rows in gamma matrix
      for (i in g_start:g_end){
        row <- strsplit(iq_file[[i]]," ")[[1]] # split the rows on the (large amount of) " "'s in the middle
        row <- row[str_detect(row,"([0-9])")] # take only the numeric elements of the vector
        g1 <- c(g1,as.numeric(row[1])) # add the values to the columns
        g2 <- c(g2,as.numeric(row[2]))
        g3 <- c(g3,as.numeric(row[3]))
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
    # Extract state frequencies
    sf1      <- strsplit(iq_file[[grep("State frequencies:",iq_file)]],":")[[1]][2]
    
    # make a list of the rows for the output dataframe
    names <- c("file_name","sequence_type","n_taxa","n_sites",num_names,"substitution_model","model_of_rate_heterogeneity","model_of_rate_heterogeneity_line2_name",
               "model_of_rate_heterogeneity_line2_value","state_frequencies")
    # Make a list of the output rows for the first output dataframe
    op <- c(op1,"amino-acid",op2,op3,num_vals,op4,mrh1,mrh2_name,mrh2,sf1)
    op_df <- data.frame(names,op, stringsAsFactors = FALSE)
    names(op_df) <- c("parameter","value")
    
    # Check whether a gamma matrix is needed
    mrh1_check <- gsub(" ","",mrh1)
    if (mrh1_check=="Uniform"){
      # If the model for rate heterogeneity is uniform, don't need to create a matrix for discrete gamma rate categories
      g_df <- "Uniform"
    } else {
      #Create the matrix for discrete gamma categories
      g_start <- grep(" Category",iq_file)+1 # get the index for the first line of the gamma categories matrix
      empty   <- which(iq_file=="") # get indexes of all empty lines
      empty   <- empty[empty>g_start] # get empty lines above gamma categories matrix
      g_end   <- empty[1]-1 # get end index for gamma categories matrix (one less than next empty line)
      end_line <- iq_file[g_end]
      # if the end isn't an empty line, subtract one from the end count 
      # to exclude lines like "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category."
      # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
      check_line <- length(strsplit(strsplit(end_line, "        " )[[1]][1]," ")[[1]])
      if (check_line > 3){
        # If the check_line is longer than 3 objects, it won't be a group for the gamma categories but an instruction
        # Eg. the relative rates line has 17 objects, because each word in the sentence is an object
        # Instructions can be excluded from the gamma matrix (but categories can't)
        g_end = g_end - 1
      }
      # initialise columns to store data in
      g1 <- c()
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
    
    # Check whether state frequencies are needed
    sf1_squashed <- gsub(" ","",sf1)
    if (sf1_squashed == "(empiricalcountsfromalignment)"){
      # Get starting line for frequencies
      start_ind <- grep("State frequencies:",iq_file) + 2
      # Take the 20 lines containing AA frequencies
      freq_lines <- iq_file[start_ind:(start_ind+19)]
      # Split up the frequency lines into the label and the frequency
      freq_split <- unlist(strsplit(freq_lines,"="))
      # Get the frequency
      freq_nums <- freq_split[c(FALSE,TRUE)]
      # Remove any spaces (from IQTree formatting)
      freq_nums <- gsub(" ","",freq_nums)
      # Get corresponding AA letter
      freq_names <- freq_split[c(TRUE,FALSE)]
      # Remove IQTree formatting
      freq_names <- gsub("pi\\(","",freq_names)
      freq_names <- gsub("\\)","",freq_names)
      freq_names <- gsub(" ","", freq_names)
      # Create a nice dataframe
      f_df <- data.frame("amino_acid" = freq_names,
                         "frequency" = freq_nums,
                         stringsAsFactors = FALSE)
    } else {
      f_df <- "State frequencies from model"
    }
    
    # Create a list of the dataframes - this will be the output
    params <- list(op_df,g_df,f_df)
    # Name the parameters so they're easy to access once you've outputted the data
    names(params) <- c("parameters","gamma_categories", "frequency")
  }
  
  # Now the information has been collected, create an output dataframe
  # make the output dataframe
  return(params)
}





phylo.collate.bootstrap <- function(alignment_folder, exec_paths, tree_folder ){
  # set the working directory to the alignment folder
  setwd(alignment_folder)
  
  # Open the original alignment results
  if (file.exists(paste0(alignment_folder,"testStatistics.csv")) == TRUE){
    alignment_df <- read.csv(paste0(alignment_folder,"testStatistics.csv"), stringsAsFactors = FALSE)
  } else if (file.exists(paste0(alignment_folder,"testStatistics.csv")) == FALSE){
    params_file <- paste0(alignment_folder,"params.csv")
    row <- read.csv(params_file,stringsAsFactors = FALSE)
    phylo.fixedtrees.run1sim(row, program_paths, tree_folder)
    alignment_df <- read.csv(paste0(alignment_folder,"testStatistics.csv"), stringsAsFactors = FALSE)
  }
  
  alignment_df$bootstrap_id <- "alignment"
  # Collect the PHI and 3Seq P-Values from the alignment df BEFORE pruning it
  PHI_sig <- alignment_df$PHI_sig[1]
  seq_sig <- alignment_df$X3SEQ_p_value[1]
  # Extract only the columns you want
  cols <- c("method", "bootstrap_id", "n_taxa", "n_sites", "tree_age", "tree1", "proportion_tree1", "tree2", "proportion_tree2", "id", "PHI_mean", "PHI_variance",
            "PHI_observed" ,"X3SEQ_num_recombinant_triplets", "X3SEQ_num_distinct_recombinant_sequences", "LM_prop_resolved_quartets","tree_proportion",
            "mean_delta_q","median_delta_q","mode_delta_q")
  alignment_df <- alignment_df[,cols]
  
  # Collate bootstrap csvs
  # Collect all the folders within the directory
  folder_paths <- list.files(alignment_folder)
  # Get the positions at which the simulations containing the id are at
  inds <- grep("bootstrap", folder_paths)
  bootstrap_paths <- folder_paths[inds] # get everything that has the word "bootstrap" in it
  bootstrap_csv_paths <- bootstrap_paths[grep("csv",bootstrap_paths)] # get all csv files
  bootstrap_folders <- setdiff(bootstrap_paths,bootstrap_csv_paths) # remove the csv files from the bootstrap folders list
  # Get the paths to the csv files
  csv_paths <- paste0(alignment_folder,bootstrap_folders,"/bootstrap_testStatistics.csv")
  # Collect all the csv files
  collate_list <- lapply(csv_paths, read.csv, stringsAsFactors = FALSE)
  # Turn the list into a dataframe
  collate_df <- Reduce(rbind, collate_list)
  collate_df <- collate_df[,cols] # extract only the cols of interest
  
  # Attach the two dataframes together
  p_value_df <- rbind(alignment_df, collate_df)
  # Output the p-value df
  bs_collated_csv <- paste0(alignment_folder,"collated_bootstrap_testStatistics.csv")
  write.csv(p_value_df,file = bs_collated_csv)
  
  # Calculate the p-values for each test statistic
  prop_resolved_quartets_sig <- calculate.p_value(p_value_df$LM_prop_resolved_quartets, p_value_df$bootstrap_id)
  nn_trimmed_sig <- calculate.p_value(p_value_df$tree_proportion, p_value_df$bootstrap_id)
  mean_delta_q_sig  <- calculate.p_value(p_value_df$mean_delta_q, p_value_df$bootstrap_id)
  median_delta_q_sig <- calculate.p_value(p_value_df$median_delta_q, p_value_df$bootstrap_id)
  mode_delta_q_sig <- calculate.p_value(p_value_df$mode_delta_q, p_value_df$bootstrap_id)
  
  # Create an output dataframe of just P-values
  op_row <- c(alignment_df[["n_taxa"]],alignment_df[["n_sites"]],alignment_df[["tree_age"]],alignment_df[["tree1"]],alignment_df[["proportion_tree1"]],alignment_df[["tree2"]],
              alignment_df[["proportion_tree2"]], alignment_df[["id"]],PHI_sig, seq_sig, prop_resolved_quartets_sig, 
              nn_trimmed_sig, mean_delta_q_sig, median_delta_q_sig, mode_delta_q_sig)
  output_df <- data.frame(matrix(nrow=0,ncol=19)) # make somewhere to store the results
  output_df <- rbind(output_df,op_row,stringsAsFactors = FALSE) # place row in dataframe
  names(output_df) <- c("n_taxa","n_sites","tree_age","tree1","proportion_tree1","tree2","proportion_tree2","id","PHI_p_value",
                        "3Seq_p_value","LM_p_value","tree_proportion_p_value","mean_delta_q_p_value",
                        "median_delta_q_p_value","mode_delta_q_p_value")
  p_value_csv <- paste0(alignment_folder,"p_value.csv")
  write.csv(output_df,file = p_value_csv)
  }





# Given two vectors (one of test statistic values, and one of ids), calculates the p-value for that alignment
calculate.p_value <- function(value_vector,id_vector){
  p_value_df <- data.frame(value_vector,id_vector, stringsAsFactors = FALSE)
  names(p_value_df) <- c("value","id")
  alignment_value <- p_value_df[which(p_value_df$id == "alignment"),1]
  if (is.na(alignment_value) == TRUE){
    # If the alignment value is NA, can't calculate a score
    # If it wasn't possible to calculate a score (will usually be a PHI score) for the alignment, output an NA
    p_value_2tail <- NA
  } else {
    # If it's possible to calculate a p-value, calculate one
    # Exclude NA rows
    p_value_df <-  p_value_df[!is.na(p_value_df$value),]
    # Find the number of bootstrap replicates and where the actual alignment value is located
    num_rows <- nrow(p_value_df) # number of bootstrap replicates + alignment value
    p_value_df <- p_value_df[order(p_value_df$value),] # order values from smallest to largest
    alignment_row <- which(p_value_df$id == "alignment") # find the ranking of the alignment value
    alignment_value <- p_value_df[alignment_row,1] # find the alignment's test statistic value
    # check whether there are other values that are the same as the alignment value
    identical_df <- subset(p_value_df,value == alignment_value)
    # if there are identical values, you don't know where the alignment actually falls within that list
    if (nrow(identical_df)>1){
      # get all the indexes of identical values
      identical_inds <- grep(alignment_value,p_value_df$value)
      # pick an ind at random
      random_identical_row <- sample(identical_inds,1)
      # For left tail probability: want to find the number of observations less than or equal to the alignment value, then divide by the number of bootstrap observations
      # If there are 8 rows and the alignment is the 5th row, then there will be 5 alignments less than or equal to the alignment value
      p_value_left <- random_identical_row/num_rows
      # For right tail probability: want to find the number of observations greater than or equal to the alignment value, then divide by the number of bootstrap observations
      # If there are 8 rows and the alignment is the 5th row, then there will be (8-5)=3 alignments greater than or equal to the alignment value
      p_value_right <- (num_rows - random_identical_row)/num_rows
    } else if (nrow(identical_df) == 1){
      # else, simply calculate the p value using the formula 
      # For left tail probability: want to find the number of observations less than or equal to the alignment value, then divide by the number of bootstrap observations
      # If there are 8 rows and the alignment is the 5th row, then there will be 5/8 alignments less than or equal to the alignment value
      p_value_left <- alignment_row/num_rows
      # For right tail probability: want to find the number of observations greater than or equal to the alignment value, then divide by the number of bootstrap observations
      # If there are 8 rows and the alignment is the 5th row, then there will be (8-5)=3 alignments greater than or equal to the alignment value
      p_value_right <- (num_rows-alignment_row)/num_rows
    }
    # To find two tailed probability, multiply the lower of those values by 2
    p_value_2tail <- 2*min(p_value_left,p_value_right) 
  }

  # return the p-value
  return(p_value_2tail)
}





# Function to select a subset of alignments from an experiment run to apply the parametric bootstrap to
# Only want to run a sample of e.g. ten alignments per point for the bootstrap so need to extract those alignments
# write a quick function to return either alignment folder (if it's one of the ten reps i.e. the last number is 10 or less) or nothing
reject.high.reps <- function(alignment_folder,max_rep = 20){
  max_rep <- as.numeric(max_rep) # make sure it's a number
  # Extract the rep number for an alignment folder
  ss <- strsplit(alignment_folder,"_")[[1]]
  rep <- ss[length(ss)] # get final object from folder name (this is where the rep number is stored in the folder name)
  rep <- gsub("/","",rep) # remove the terminal slash from the rep number
  rep <- as.numeric(rep) # convert rep number to numeric
  # If the rep is one of the ones you want, return the name for collection
  if (is.na(rep) == TRUE || (rep =="") == TRUE){
    # if rep is NA, return NULL
    return(NULL)
  } else if (rep <= max_rep){
    # if rep in 1:max_rep, return folder to run bootstraps
    return(alignment_folder)
  } else if (rep > max_rep){
    # if rep > max_rep, ignore folder
    return(NULL)
  }
}





# Function to select a subset of alignments from an experiment run to apply the parametric bootstrap to
# Only want to run the alignment if the proportion of introgressed DNA is 0%,10%,20%,30%,40% or 50%
# write a quick function to return either alignment folder (if the proportion of introgressed DNA is a multiple of 10%) or nothing
reject.exp3.bootstrap.run <- function(alignment_folder){
  # split the alignment folder by the underscores
  split_f <- strsplit(alignment_folder, "_")[[1]]
  # find where the first part of the name is (look for the "FixedTrees" term) and get that bit onwards
  split_f <- split_f[grep("FixedTrees",split_f):length(split_f)]
  # the proportion of tree 2 will be after the tree 2 structure
  K <- as.numeric(split_f[6])
  valid_K <- c(0,0.1,0.2,0.3,0.4,0.5)
  # check if this alignment needs a parametric bootstrap
  if (isTRUE(all.equal(K,0))||isTRUE(all.equal(K,0.1))||isTRUE(all.equal(K,0.2))||isTRUE(all.equal(K,0.3))||isTRUE(all.equal(K,0.4))||isTRUE(all.equal(K,0.5))){
    # if it does, return it 
    return(alignment_folder)
  }
}




