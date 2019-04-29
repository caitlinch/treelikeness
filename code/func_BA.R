# Functions to run test statistics on empirical datasets

library(phytools)
library(ape)
library(phangorn)

empirical.runTS <- function(alignment_path, program_paths, bootstrap_id = FALSE){
  # extract the alignment folder from the alignment path
  alignment_folder <- paste0(dirname(alignment_path),"/")
  # Extract the dataset name (basename of alignment folder: element after last "/" in alignment_folder)
  dataset <- basename(alignment_folder)
  # Create some folder and filenames
  loci_name <- gsub(".nex","",basename(alignment_path))
  log_folder <- paste0(alignment_folder,loci_name,"/")
  
  # If the log file doesn't exist, create it 
  if (dir.exists(log_folder) == FALSE){
    dir.create(log_folder) # create a new folder to store the log files from the executables for this loci in
  }
  
  # Open the nexus file and get the number of taxa and the number of characters 
  n <- read.nexus.data(alignment_path)
  n_taxa <- length(n)
  n_char <- length(unlist(n[1]))
  
  # Run IQ-tree on the alignment (if it hasn't already been run), and get the likelihood mapping results
  call.IQTREE.quartet(program_paths[["IQTree"]],alignment_path,n_taxa)
  
  # Change to the folder for this alignment - means that 3seq and Phi files will be saved into a unique folder
  setwd(log_folder)
  # Get paths to PhiPac, 3SEQ
  phi_path <- program_paths[["Phi"]] # get path to phipack executable
  seq_path <- program_paths[["3seq"]] # get path to 3seq executable
  # Phipack only reads in Phylip or fasta format - need to convert if the alignment is a nexus file (using the nexus data opened above)
  fasta.name <- paste0(log_folder,loci_name,".fasta") # make a name for the fasta alignment by adding .fasta (super original ;) )
  write.fasta(sequences = n,names = names(n), file.out = fasta.name) # output alignment as a fasta format
  # run PHIPACK and 3seq 
  # Note that Phi and 3Seq will only be run if they haven't already been run (checks for log files)
  if (file.exists(paste0(log_folder,"Phi.log")) == FALSE){
    phi_command <- paste0(phi_path," -f ",fasta.name, " -v") # assemble system command
    system(phi_command) #call phipack
  }
  if (file.exists(paste0(log_folder,"3s.log")) == FALSE){
    seq_command <- paste0(seq_path," -f ", fasta.name)
    system(seq_command) #call 3SEQ
  }
  
  # Both PhiPack and 3SEQ will have created a file with the information about the statistics in
  # Extract significance from Phi Pack output
  phi_file <- paste0(log_folder,"Phi.log")
  phi_file <- readLines(phi_file)
  ind      <- grep("p-Value",phi_file)
  phi_sig <- as.numeric(strsplit(phi_file[ind+3],":")[[1]][2])
  ind      <- grep("PHI Values",phi_file)
  phi_mean <- as.numeric(strsplit(phi_file[ind+4],"      ","")[[1]][2])
  phi_var <- as.numeric(strsplit(phi_file[ind+5],"      ","")[[1]][2])
  phi_obs <- as.numeric(strsplit(phi_file[ind+6],"      ","")[[1]][2])
  
  # Extract results output from 3Seq output
  seq_file <- paste0(log_folder,"3s.log")
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
  
  # Change back to directory containing alignments and iqtree files
  setwd(alignment_folder)
  
  # Extract quartet mapping
  iq_log_path <- paste0(alignment_path,".iqtree")
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
  pdmm <- as.matrix(mldist.pdm(alignment_path)) # open pairwise distance matrix as a matrix
  deltaplot_results <- delta.plot(pdmm, k = 51, plot = FALSE) # calculate the delta.plot
  counts <- deltaplot_results$counts
  intervals <- seq(0,1,(1/(length(counts)-1)))
  deltaplot_df <- data.frame(intervals,counts)
  names(deltaplot_df) <- c("intervals","counts")
  deltaplot_df_name <- paste0(log_folder,loci_name,"_deltaplot_histogram.csv")
  write.csv(deltaplot_df,file = deltaplot_df_name)
  # Want to calculate the mean and median delta q value - unfortunately the delta.plot function doesn't output raw data, so make a pseudo data set using the histogram values
  mean_dq <- mean(rep(deltaplot_df$intervals,deltaplot_df$counts)) # turn the interval data into a long list of "raw" values and calculate the mean
  median_dq <- median(rep(deltaplot_df$intervals,deltaplot_df$counts)) # turn the interval data into a long list of "raw" values and calculate the median
  mode_dq <- deltaplot_df[order(deltaplot_df$count, decreasing = TRUE),][1,1] # sort the dataframe by count values and extract the mode
  
  # Run my test statistics
  # run pdm ratio (TS1) (modified splittable percentage)
  splittable_percentage <- pdm.ratio(iqpath = program_paths[["IQTree"]], path = alignment_path)
  # run normalised.pdm.difference.sum (TS2a) (sum of difference of normalised matrix)
  npds <- normalised.pdm.diff.sum(iqpath = program_paths[["IQTree"]], path = alignment_path)
  # run normalised pdm difference average (TS2b) (mean of difference of normalised matrix)
  npdm <- normalised.pdm.diff.mean(iqpath = program_paths[["IQTree"]], path = alignment_path)
  
  # Run trimmed and untrimmed versions of the split decomposition and NeighborNet tree proportion
  # Create a new nexus file with a taxa block
  new_nexus_file <- paste0(alignment_folder,loci_name,"_withTaxaBlock.nexus")
  write.nexus.data(n, file = new_nexus_file,datablock = FALSE, interleaved = FALSE)
  # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
  nexus_edit <- readLines(new_nexus_file) # open the new nexus file
  ind <- grep("BEGIN CHARACTERS",nexus_edit)+2 # find which line
  nexus_edit[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
  writeLines(nexus_edit,new_nexus_file) # output the edited nexus file
  # Call the test statistic functions
  initial_iqtree_tree <- paste0(alignment_path,".treefile")
  sd_untrimmed <- tree.proportion(iqpath = program_paths[["IQTree"]], splitstree_path = program_paths[["SplitsTree"]], path = new_nexus_file, network_algorithm = "split decomposition", trimmed = FALSE, tree_path = initial_iqtree_tree, run_IQTREE = FALSE)
  nn_untrimmed <- tree.proportion(iqpath = program_paths[["IQTree"]], splitstree_path = program_paths[["SplitsTree"]], path = new_nexus_file, network_algorithm = "neighbournet", trimmed = FALSE, tree_path = initial_iqtree_tree, run_IQTREE = FALSE)
  sd_trimmed <- tree.proportion(iqpath = program_paths[["IQTree"]], splitstree_path = program_paths[["SplitsTree"]], path = new_nexus_file, network_algorithm = "split decomposition", trimmed = TRUE, tree_path = initial_iqtree_tree, run_IQTREE = FALSE)
  nn_trimmed <- tree.proportion(iqpath = program_paths[["IQTree"]], splitstree_path = program_paths[["SplitsTree"]], path = new_nexus_file, network_algorithm = "neighbournet", trimmed = TRUE, tree_path = initial_iqtree_tree, run_IQTREE = FALSE)
  
  # Name the test statistics file (if it's a  bootstrap replicate, add the replicate number!)
  # Create the variable for bootstrap replicate
  if (bootstrap_id == FALSE){
    results_file <- paste0(alignment_folder,loci_name,"_testStatistics.csv")
    rep_id = "empiricalAlignment"
  } else {
    results_file <- paste0(alignment_folder,loci_name,"_bootstrapReplicate",bootstrap_id,"_testStatistics.csv")
    rep_id <- paste0("bootstrapReplicate",bootstrap_id)
  }
  
  # Make somewhere to store the results
  df_names <- c("dataset","loci","bootstrap_replicate_id","n_taxa","n_sites","alignment_file",
                "PHI_mean","PHI_variance","PHI_observed","PHI_sig","3SEQ_num_recombinant_triplets","3SEQ_num_distinct_recombinant_sequences","3SEQ_p_value","num_quartets",
                "num_resolved_quartets","prop_resolved_quartets","num_partially_resolved_quartets","num_unresolved_quartets", "splittable_percentage","pdm_difference",
                "pdm_average","split_decomposition_untrimmed", "neighbour_net_untrimmed", "split_decomposition_trimmed","neighbour_net_trimmed","mean_delta_q","median_delta_q","mode_delta_q")
  df <- data.frame(matrix(nrow=0,ncol=length(df_names))) # create an empty dataframe of the correct size
  op_row <- c(dataset,loci_name,rep_id,n_taxa,n_char,alignment_path,
              phi_mean,phi_var,phi_obs,phi_sig,num_trips,num_dis,seq_sig,total_q,resolved_q,
              prop_resolved,partly_resolved_q,unresolved_q,splittable_percentage,npds,npdm,sd_untrimmed,nn_untrimmed,
              sd_trimmed,nn_trimmed,mean_dq,median_dq,mode_dq) # collect all the information
  df <- rbind(df,op_row,stringsAsFactors = FALSE) # place row in dataframe
  names(df) <- df_names # add names to the df so you know what's what
  
  write.csv(df,file = results_file)
}

empirical.parametric.bootstrap <- function(){
  # If it hasn't already been run, call and run IQTree
  call.IQTREE(program_paths["IQTree"],alignment_path)
  
  #Extract the parameters from the .iqtree log file.
  params <- get.simulation.parameters(paste0(alignment_path,".iqtree"))
  
  # output_dataframe: dataset, loci name, number of taxa, number of characters, test statistic values
}