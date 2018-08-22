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
phylo.parametric.bootstrap <- function(alignment_folder,n_reps,iq_path,splitstree_path, phipack_path, threeseq_path) {
  # extract the name of the .iqtree file that contains the parameters for the simulation
  dotiqtree_path <- paste0(alignment_folder, list.files(alignment_folder)[grep("iqtree", list.files(alignment_folder))])
  params <- get.simulation.parameters(dotiqtree_path) #need to feed in .iqtree file
  
  # Open the ML tree from IQ-Tree
  # extract the name of the file that contains the ML tree calculated by IQ-Tree for the original alignment (the alignment with recombination)
  all_tree_paths <- paste0(alignment_folder, list.files(alignment_folder)[grep("treefile", list.files(alignment_folder))]) #get all tree files (will be three - one from IQ-Tree, and one for each tree 1 and tree 2)
  ML_tree_path <- all_tree_paths[grep("alignment", all_tree_paths)] # get only the treefile for the alignment
  ML_tree <- read.tree(ML_tree_path)
  
  # fill out the rep numbers (padded with 0s to get to 4 digits)
  rep_ids <- 1:n_reps
  rep_ids <- sprintf("%04d", rep_ids)
  lapply(rep_ids,do.1.bootstrap,params,ML_tree,alignment_folder,iq_path,splitstree_path, phipack_path, threeseq_path)
}

# Given the relevant information, run one parametric bootstrap (create the alignment and run the test statistics, output the p-values as a vector)
do.1.bootstrap <- function(rep_number,params,tree,alignment_folder,iq_path,splitstree_path, phipack_path, threeseq_path){
  # Create a new folder name to store this alignment and its outputs in
  bs_folder <- paste0(alignment_folder,"bootstrap_",rep_number,"/")
  
  if (dir.exists(bs_folder)==TRUE){
    redo <- TRUE
  } else if (dir.exists(bs_folder)==FALSE){
    # if the folder doesn't exist, create it
    dir.create(bs_folder)
    redo <- FALSE  # if the directory and the alignment don't exist, this is not a redo and the alignment needs to be tested
  } 
  
  # Only create alignments, run IQ-Tree and all the tests/test statistics and everything if this is NOT a redo (redo == FALSE)
  if (redo == FALSE){
    # Change the working directory to the bootstrap folder
    setwd(bs_folder)
    
    # Create the alignment 
    # Sample code for generating a parametric DNA sequence if you have a tree
    # s1 = simSeq(t1, l = 500, type="DNA", bf=c(.25,.25,.25,.25), Q=c(1,1,1,1,1,1), rate=1)
    # s2 = simSeq(t2, l = 500, type="DNA", bf=c(.25,.25,.25,.25), Q=c(1,1,1,1,1,1), rate=1)
    # aln = c(s1, s2)
    
    # Extract the parameters you need to enter into simSeq
    n_bp = as.numeric(params$parameters[4,2]) # want a sequence that is 1300 base pairs long
    # Extract the vector form of the rate matrix 
    m <- params$Q_rate_matrix[,2:5] # extract the square block with the rates and not the header column
    Q_vec <- c(m[2,1],m[3,1],m[3,2],m[4,1],m[4,2],m[4,3]) # extract the rates 
    # Extract the base frequencies in the following order: A, C, G, T
    base_freqs <- c(as.numeric(params$parameters[[12,2]]), as.numeric(params$parameters[[13,2]]), as.numeric(params$parameters[[14,2]]), as.numeric(params$parameters[[15,2]]))
    seq_type <- "DNA" # generate DNA sequence
    
    # Generate the DNA sequence
    gamma_categories <- params$gamma_categories
    if (class(gamma_categories) == "character"){
      # This means that gamma_categories = "Uniform"
      # No rate heterogeneity exists - generate the DNA all in one chunk
      dna_sim <- simSeq(tree, l = n_bp, type = seq_type, bf = base_freqs, Q = Q_vec)
    } else if (class(gamma_categories) == "data.frame") {
      # If the gamma categories are a dataframe, there's rate heterogeneity happening
      # If there's rate heterogeneity, stick together a bunch of alignments to make one big alignment
      gamma_categories$n_bp <- (gamma_categories$proportion * n_bp)
      
      # Make sure that the sum of alignment types will be n_bp (1300 bp)
      tot <- sum(gamma_categories$n_bp)
      diff <- n_bp - tot
      if (diff > 0) {
        # if there are less bp in the gamma categories than the desired number of base pairs
        # need to add the diff onto one category
        # Randomly pick one row to add
        change_row <- sample(1:nrow(gamma_categories),1)
        # Add the difference to the randomly chosen row
        gamma_categories$n_bp[change_row] <- gamma_categories$n_bp[change_row] + diff
      } else if (diff < 0){
        # if there are more bp in the gamma categories than the desired number of base pairs
        # need to subtract the diff from one category
        # Randomly pick one row to subtract
        change_row <- sample(1:nrow(gamma_categories),1)
        # Subtract the difference from the randomly chosen row
        gamma_categories$n_bp[change_row] <- gamma_categories$n_bp[change_row] - diff
      }
    # Now you have the number of bp for each rate, you can create the alignments and stick them together
    for (i in 1:nrow(gamma_categories)){
      # Iterate through each row of the gamma categories matrix
      # Extract the rate for that chunk of the sequence
      chunk_rate <- gamma_categories$relative_rate[[i]]
      # Extract the number of base pairs for that chunk of the sequence
      chunk_length <- gamma_categories$n_bp[[i]]
      # Simulate the DNA
      dna_chunk <- simSeq(tree, l = chunk_length, type = seq_type, bf = base_freqs, Q = Q_vec, rate = chunk_rate)
      dnabin_chunk <- as.DNAbin(dna_chunk)
      # Stick the DNA together with the DNA already created (or the empty object to store DNA in)
      if (exists("dna_bin") == FALSE){
        # If the dna_bin object is empty, this is the first chunk. Make it the dna_bin object
        dna_bin <- dnabin_chunk
      } else if (exists("dna_bin") == TRUE){
        # Otherwise, the dna_bin object exists. Concatenate the dna chunk onto the dna_bin object
        dna_bin <- cbind(dna_bin,dnabin_chunk, check.names = TRUE, fill.with.gaps = TRUE, quiet = FALSE) # concatenate the two alignments
      }
    }
    dna_sim <- as.phyDat(dna_bin) # convert the alignment to  phydat format
    }
      
    # Save the DNA alignment
    bs_al <- paste0(bs_folder, "alignment.nexus")
    write.phyDat(dna_sim,file = bs_al, format = "nexus",interleaved = TRUE, datablock = FALSE) # write the output as a nexus file
    # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
    nexus <- readLines(bs_al) # open the new nexus file
    ind <- grep("BEGIN CHARACTERS",nexus)+2 # find which line
    nexus[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
    writeLines(nexus,bs_al) # output the edited nexus file
    
    # The alignment now definitely exists. Now you can run IQ-tree on the alignment
    n_taxa <- as.numeric(params$parameters[3,2]) # extract the number of taxa from the parameters 
    call.IQTREE.quartet(iq_path,bs_al,n_taxa)
  }
  bs_al <- paste0(bs_folder, "alignment.nexus") # alignment name
  
  ## Calculate the test statistics
  # run PHIPACK and 3seq
  phi_path <- phipack_path # get path to phipack executable
  seq_path <- threeseq_path # get path to 3seq executable
  filetype = tail(strsplit(bs_al,"\\.")[[1]],n=1) # extract file format
  if (filetype == "fasta"){
    # if the alignment is already in fasta format, run PhiPack through R
    phi_command <- paste0(phi_path," -f ",bs_al, " -v") # assemble system command
    system(phi_command) #call phipack
    
    seq_command <- paste0(seq_path," -f ", bs_al)
    system(seq_command) #call 3SEQ
  } else if (filetype == "nexus"){
    # Phipack only reads in Phylip or fasta format - need to convert if the alignment is a nexus file
    data = read.nexus.data(bs_al) # read in nexus format alignment
    fasta.name <- paste0(bs_al,".fasta") # make a name for the fasta alignment by adding .fasta (super original ;) )
    write.fasta(sequences = data,names = names(data), file.out = fasta.name) # output alignment as a fasta format
    phi_command <- paste0(phi_path," -f ",fasta.name, " -v") # assemble system command as above
    system(phi_command) # run PHI test on the new fasta alignment
    
    seq_command <- paste0(seq_path," -f ", fasta.name)
    system(seq_command) #call 3SEQ
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
  
  # Run my test statistics
  # Want the path to be the path to the alignment you're testing - here that's the bootstrap alignment (bs_al)
  # run pdm ratio (TS1) (modified splittable percentage)
  splittable_percentage <- pdm.ratio(iqpath = iq_path, path = bs_al)
  # run normalised.pdm.difference.sum (TS2a) (sum of difference of normalised matrix)
  npds <- normalised.pdm.diff.sum(iqpath = iq_path, path = bs_al)
  # run normalised pdm difference average (TS2b) (mean of difference of normalised matrix)
  npdm <- normalised.pdm.diff.mean(iqpath = iq_path, path = bs_al)
  # run split decomposition (TS3) (split decomposition using splitstree)
  sd <- SplitsTree.decomposition.statistic(iqpath = iq_path, splitstree_path = splitstree_path, path = bs_al, network_algorithm = "split decomposition")
  # run NeighbourNet (TS3, with neighbour net not split decomposition using splitstree)
  nn <- SplitsTree.decomposition.statistic(iqpath = iq_path, splitstree_path = splitstree_path, path = bs_al, network_algorithm = "neighbournet")
  
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
  params_csv$num_quartets <- total_q
  params_csv$num_resolved_quartets <- resolved_q
  params_csv$prop_resolved_quartets <- prop_resolved
  params_csv$num_partially_resolved_quartets <- partly_resolved_q
  params_csv$num_unresolved_quartets <- unresolved_q
  params_csv$splittable_percentage <- splittable_percentage
  params_csv$pdm_difference <- npds
  params_csv$pdm_average <- npdm
  params_csv$split_decomposition <- sd
  params_csv$neighbour_net <- nn
  params_csv$bootstrap_id <- paste0("bootstrap_",rep_number)
  
  # Output results as a csv in the bootstrap folder
  bs_results_csv <- paste0(bs_folder,"bootstrap_testStatistics.csv")
  write.csv(params_csv,file = bs_results_csv)
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
      # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
      check_line <- gsub(" ","",(strsplit(end_line, "        " )[[1]][1]))
      if (length(check_line) > 3){
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
    # Name the parameters so they're easy to access once you've outputted the data
    names(params) <- c("parameters","gamma_categories")
  }
  
  # Now the information has been collected, create an output dataframe
  # make the output dataframe
  return(params)
}

phylo.collate.bootstrap <- function(alignment_folder){
  # set the working directory to the alignment folder
  setwd(alignment_folder)
  
  # Open the original alignment results
  alignment_df <- read.csv(paste0(alignment_folder,"testStatistics.csv"))
  alignment_df$bootstrap_id <- "alignment"
  # Collect the PHI and 3Seq P-Values from the alignment df BEFORE pruning it
  PHI_sig <- alignment_df$PHI_sig[1]
  seq_sig <- alignment_df$X3SEQ_p_value[1]
  # Extract only the columns you want
  cols <- c("method", "bootstrap_id", "n_taxa", "n_sites", "tree_age", "tree1", "proportion_tree1", "tree2", "proportion_tree2", "id", "PHI_mean", "PHI_variance",
            "PHI_observed" ,"X3SEQ_num_recombinant_triplets", "X3SEQ_num_distinct_recombinant_sequences", "prop_resolved_quartets", "splittable_percentage",
            "pdm_difference", "pdm_average", "split_decomposition", "neighbour_net")
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
  PHI_observed_sig <- calculate.p_value(p_value_df$PHI_mean, p_value_df$bootstrap_id)
  PHI_observed_sig <- calculate.p_value(p_value_df$PHI_observed, p_value_df$bootstrap_id)
  x3seq_sig <- calculate.p_value(p_value_df$X3SEQ_num_distinct_recombinant_sequences, p_value_df$bootstrap_id)
  prop_resolved_quartets_sig <- calculate.p_value(p_value_df$prop_resolved_quartets, p_value_df$bootstrap_id)
  splittable_percentage_sig <- calculate.p_value(p_value_df$splittable_percentage, p_value_df$bootstrap_id)
  pdm_difference_sig <- calculate.p_value(p_value_df$pdm_difference, p_value_df$bootstrap_id)
  pdm_average_sig <- calculate.p_value(p_value_df$pdm_average, p_value_df$bootstrap_id)
  split_decomposition_sig <- calculate.p_value(p_value_df$split_decomposition, p_value_df$bootstrap_id)
  neighbour_net_sig <- calculate.p_value(p_value_df$neighbour_net, p_value_df$bootstrap_id)
  
  # Create an output dataframe of just P-values
  op_row <- c(alignment_df[["n_taxa"]],alignment_df[["n_sites"]],alignment_df[["tree_age"]],alignment_df[["tree1"]],alignment_df[["proportion_tree1"]],alignment_df[["tree2"]],
              alignment_df[["proportion_tree2"]], alignment_df[["id"]],PHI_sig, PHI_mean_sig, PHI_observed_sig, seq_sig, x3seq_sig, prop_resolved_quartets_sig, splittable_percentage_sig, pdm_difference_sig, pdm_average_sig, 
              split_decomposition_sig, neighbour_net_sig)
  output_df <- data.frame(matrix(nrow=0,ncol=19)) # make somewhere to store the results
  output_df <- rbind(output_df,op_row,stringsAsFactors = FALSE) # place row in dataframe
  names(output_df) <- c("n_taxa","n_sites","tree_age","tree1","proportion_tree1","tree2","proportion_tree2","id","PHI_p_value","PHI_mean_p_value","PHI_observed_p_value",
                        "3Seq_p_value","num_recombinant_sequences_p_value","likelihood_mapping_p_value","splittable_percentage_p_value","pdm_difference_p_value",
                        "pdm_average_p_value","split_decomposition_p_value","neighbour_net_p_value")
  p_value_csv <- paste0(alignment_folder,"p_value.csv")
  write.csv(output_df,file = p_value_csv)
}

# Given two vectors (one of test statistic values, and one of ids), calculates the p-value for that alignment
calculate.p_value <- function(value_vector,id_vector){
  p_value_df <- data.frame(value_vector,id_vector, stringsAsFactors = FALSE)
  names(p_value_df) <- c("value","id")
  # Order by test statistic value
  p_value_df <- p_value_df[order(p_value_df$value),]
  # Find the number of bootstrap replicates and where the actual alignment value is located
  num_bs <- nrow(p_value_df)-1 # number of bootstrap replicates - need to subtract 1 because this includes the alignment value 
  alignment_row <- which(p_value_df$id == "alignment") # find the ranking of the alignment value
  # For left tail probability: want to find the number of observations less than or equal to the alignment value, then divide by the number of bootstrap observations
  n_obs_smaller <- alignment_row - 1
  p_value_left <- n_obs_smaller/num_bs
  # For right tail probability: want to find the number of observations greater than or equal to the alignment value, then divide by the number of bootstrap observations
  n_obs_larger <- num_bs - n_obs_smaller
  p_value_right <- n_obs_larger/num_bs
  # To find two tailed probability, multiply the lower of those values by 2
  p_value_2tail <- 2*min(p_value_left,p_value_right) 
  # return the p-value
  return(p_value_2tail)
}


