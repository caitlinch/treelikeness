# Functions to run test statistics on empirical datasets

library(phytools)
library(ape)
library(phangorn)

empirical.runTS <- function(alignment_path, program_paths, bootstrap_id){
  print("in empirical.runTS")
  print(alignment_path)
  # extract the alignment folder from the alignment path
  alignment_folder <- paste0(dirname(alignment_path),"/")
  output_id <- gsub(".nex","",basename(alignment_path))
  # Create some folder and filenames
  if (bootstrap_id == "alignment"){
    # Get the alignment name and remove the extension to get the loci name
    loci_name <- gsub(".nex","",basename(alignment_path))
    # Extract the dataset name (basename of alignment folder: element after last "/" in alignment_folder)
    dataset <- basename(alignment_folder)
  } else {
    # If the alignment is a bootstrap replicate, need to remove the bootstrap rep number to get the loci name
    loci_name <- gsub(".nex","",basename(alignment_path)) # get the basis of the loci name
    loci_list <- unlist(strsplit(loci_name, "_")) # break the alignment name into bits
    max_ind <- grep("bootstrapReplicate",loci_list) - 1 # which ind is the bootstrapReplicate at?
    loci_list <- loci_list[1:max_ind] # get only parts of alignment name
    loci_name <- paste(loci_list,collapse="_") # squash the alignment name together
    # Extract the dataset name (basename of alignment folder: element after second-last "/" in alignment_folder)
    dataset <- basename(dirname(alignment_folder))
  }
  
  if (bootstrap_id == "alignment"){
    # If this is not a bootstrap replicate, you need to create a folder to store the program logs in
    # Otherwise they will get overwritten for each locus you run - this means you can keep them for late
    log_folder <- paste0(alignment_folder,loci_name,"/")
    # make a rep id - e.g. an id to store in output df so you can identify what's a simulation and what's from empirical data
    rep_id <- loci_name
  } else {
    # If the run is a bootstrap replicate, you can just save the information in its folder (only 1 alignment per folder in bootstrap replicate folders)
    log_folder <- alignment_folder
    # make a rep id - e.g. an id to store in output df so you can identify what's a simulation and what's from empirical data
    rep_id <- bootstrap_id
  }
  
  
  # If the log file doesn't exist, create it 
  if (dir.exists(log_folder) == FALSE){
    dir.create(log_folder) # create a new folder to store the log files from the executables for this loci in
  }
  
  # Open the nexus file and get the number of taxa and the number of characters 
  n <- read.nexus.data(alignment_path)
  n_taxa <- length(n)
  n_char <- length(unlist(n[1]))
  
  # Run IQ-tree on the alignment (if it hasn't already been run), and get the likelihood mapping results
  # Check that the original alignment ran ok
  if (file.exists(paste0(alignment_path,".iqtree")) == FALSE || file.exists(paste0(alignment_path,".treefile")) == FALSE || file.exists(paste0(alignment_path,".lmap.eps")) == FALSE){
    print("need to run IQ-Tree - missing files")
    n <- read.nexus.data(alignment_path)
    n_taxa <- length(n)
    # Run IQ-tree on the alignment (if it hasn't already been run), and get the likelihood mapping results
    print("run IQTree")
    call.IQTREE.quartet(program_paths[["IQTree"]],alignment_path,n_taxa)
  }
  
  # Change to the log (storage for log files) folder for this alignment - means that 3seq and Phi files will be saved into a unique folder
  print("run PHI and 3SEQ")
  setwd(log_folder)
  # Get paths to PhiPac, 3SEQ
  phi_path <- program_paths[["Phi"]] # get path to phipack executable
  seq_path <- program_paths[["3seq"]] # get path to 3seq executable
  # Phipack only reads in Phylip or fasta format - need to convert if the alignment is a nexus file (using the nexus data opened above)
  fasta.name <- paste0(log_folder,output_id,".fasta") # make a name for the fasta alignment by adding .fasta (super original ;) )
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
  print("extract likelihood mapping")
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
  print("estimate delta plots")
  pdmm <- as.matrix(mldist.pdm(alignment_path)) # open pairwise distance matrix as a matrix
  deltaplot_results <- delta.plot(pdmm, k = 51, plot = FALSE) # calculate the delta.plot
  counts <- deltaplot_results$counts
  intervals <- seq(0,1,(1/(length(counts)-1)))
  deltaplot_df <- data.frame(intervals,counts)
  names(deltaplot_df) <- c("intervals","counts")
  deltaplot_df_name <- paste0(alignment_folder,output_id,"_deltaplot_histogram.csv")
  write.csv(deltaplot_df,file = deltaplot_df_name, row.names = FALSE)
  # Want to calculate the mean and median delta q value - unfortunately the delta.plot function doesn't output raw data, so make a pseudo data set using the histogram values
  mean_dq <- mean(rep(deltaplot_df$intervals,deltaplot_df$counts)) # turn the interval data into a long list of "raw" values and calculate the mean
  median_dq <- median(rep(deltaplot_df$intervals,deltaplot_df$counts)) # turn the interval data into a long list of "raw" values and calculate the median
  mode_dq <- deltaplot_df[order(deltaplot_df$count, decreasing = TRUE),][1,1] # sort the dataframe by count values and extract the mode
  
  # Run my test statistics
  print("run my test statistics")
  # run pdm ratio (TS1) (modified splittable percentage)
  splittable_percentage <- pdm.ratio(iqpath = program_paths[["IQTree"]], path = alignment_path)
  # run normalised.pdm.difference.sum (TS2a) (sum of difference of normalised matrix)
  npds <- normalised.pdm.diff.sum(iqpath = program_paths[["IQTree"]], path = alignment_path)
  # run normalised pdm difference average (TS2b) (mean of difference of normalised matrix)
  npdm <- normalised.pdm.diff.mean(iqpath = program_paths[["IQTree"]], path = alignment_path)
  
  # Run trimmed and untrimmed versions of the split decomposition and NeighborNet tree proportion
  # Splitstree needs a specific file format - so create a new nexus file with a taxa block
  new_nexus_file <- paste0(alignment_folder,output_id,"_withTaxaBlock.nexus")
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
  
  # Name the test statistics file using the output id (this way if it's a  bootstrap replicate, it adds the replicate number!)
  print(paste0("output results for ",output_id))
  results_file <- paste0(alignment_folder,output_id,"_testStatistics.csv")
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
  
  write.csv(df,file = results_file, row.names = FALSE)
}



do1.empirical.parametric.bootstrap <- function(bootstrap_id, empirical_alignment_path, alignment_params, program_paths){
  print("in do1.empirical.parametric.bootstrap")
  print(empirical_alignment_path)
  print(bootstrap_id)
  # Create the folder for this replicate, gather and create filenames
  loci_name <- gsub(".nex","",basename(empirical_alignment_path))
  bootstrap_name <- paste0(loci_name,"_",bootstrap_id) # this will be the name of the alignment
  bootstrap_folder <- paste0(dirname(empirical_alignment_path),"/",bootstrap_name,"/") # folder to store results from this bootstrap in
  # If the bootstrap folder doesn't exist, create it
  if (dir.exists(bootstrap_folder) == FALSE){
    dir.create(bootstrap_folder)
  }
  shuffled_alignment_path <- paste0(bootstrap_folder,bootstrap_name,"_shuffled_noGaps.nex")
  bootstrap_alignment_path <- paste0(bootstrap_folder,bootstrap_name,".nex")
  empirical_alignment_tree_path <- paste0(empirical_alignment_path,".treefile")
  empirical_alignment_tree <- read.tree(empirical_alignment_tree_path)
  
  # open the empirical alignment to get information about the sequence
  print("open nexus files")
  n <- read.nexus.data(empirical_alignment_path)
  p <- phyDat(n)
  new_aln <- p
  
  # First generate completely new DNA using the params
  # Create an alignment for this replicate using the alignment params - name will be loci_bootstrapReplicateXXXX
  if (file.exists(bootstrap_alignment_path) == FALSE) {
    print("creating alignment")
    # Sample code for generating a parametric DNA sequence if you have a tree
    # s1 = simSeq(t1, l = 500, type="DNA", bf=c(.25,.25,.25,.25), Q=c(1,1,1,1,1,1), rate=1)
    # s2 = simSeq(t2, l = 500, type="DNA", bf=c(.25,.25,.25,.25), Q=c(1,1,1,1,1,1), rate=1)
    # aln = c(s1, s2) # concatenate the alignments
    
    # Extract the parameters you need to enter into simSeq
    n_bp = as.numeric(alignment_params$parameters[4,2])  # sequence length should be the same as in the original sequence
    # Extract the vector form of the rate matrix 
    m <- alignment_params$Q_rate_matrix[,2:5] # extract the square block with the rates and not the header column
    Q_vec <- c(m[2,1],m[3,1],m[3,2],m[4,1],m[4,2],m[4,3]) # extract the rates 
    # Extract the base frequencies in the following order: A, C, G, T
    base_freqs <- c(as.numeric(alignment_params$parameters[[12,2]]), as.numeric(alignment_params$parameters[[13,2]]),
                    as.numeric(alignment_params$parameters[[14,2]]), as.numeric(alignment_params$parameters[[15,2]]))
    seq_type <- "DNA" # generate DNA sequence
    
    # need to generate the number of sites for each gamma category:
    g_cat <- alignment_params$gamma_categories
    # If uniform gamma categories, can just simulate a sequence straight off
    if (class(g_cat) == "character"){
      new_aln <- simSeq(x = empirical_alignment_tree, l = n_bp, type = seq_type, bf = base_freqs, Q = Q_vec)
    } else if (class(g_cat) == "data.frame") {
      # If there are gamma categories, create multiple alignments and concatenate them: one for each gamma category
      cat_sites <- round(n_bp*g_cat$proportion,digits = 0)
      g_cat$cat_sites <- cat_sites
      g_cat <- fix.gammaCategory.siteNums(g_cat,n_bp) # make sure the sum of the gamma sites category is correct
      # Get the number of gamma rate categories
      num_g_cat <- length(g_cat$category)
      # Initialise an empty sequence
      aln_exists <- FALSE
      
      # Iterate through the gamma categories and create an alignment with each of the relative rates
      for (i in 1:num_g_cat){
        row <- g_cat[i,]
        # Generate the DNA sequence using the informtion from the parameter list and the rate and number of sites for THIS GAMMA CATEGORY
        dna_sim <- simSeq(x = empirical_alignment_tree, l = row$cat_sites, type = seq_type, bf = base_freqs, Q = Q_vec, rate = row$relative_rate)
        #concatenate the sequence with the new alignment
        if (aln_exists == FALSE){
          # If the alignment doesn't exist, create it
          new_aln <- dna_sim
          aln_exists <- TRUE
        } else if (aln_exists == TRUE){
          # If the alignment does exist, concatenate it
          new_aln <- c(new_aln,dna_sim)
        }
      }
    }
    
    # Second, randomise sites in new_aln (otherwise masking will mean the gamma categories disproportionately get affected)
    print("shuffling alignment")
    #   Turn the phyDat into a dataframe and shuffle the rows using sample()
    #   Rows are shuffled as columns represent sequences - shuffling the rows keeps the relationships between species (shown by rows) 
    #   but rearranges the order the relationships occur in
    new_aln_df <- as.data.frame(new_aln) # turn the new alignment into a dataframe
    new_aln_df <- new_aln_df[sample(nrow(new_aln_df)),] # sample the rows randomly (this will randomly distribute the gamma categories throughout)
    rownames(new_aln_df) <- 1:nrow(new_aln_df) # reset the row names as 1:nrows (they got mixed up when the sampling occurred)
    # Turn the dataframe back into an alignment:
    new_aln_shuffled <- as.phyDat(new_aln_df)
    #write the shuffled alignment to disk
    write.phyDat(new_aln_shuffled, file = shuffled_alignment_path, format = "nexus", interleave = FALSE)
    # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
    nexus_edit <- readLines(shuffled_alignment_path) # open the new nexus file
    ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
    nexus_edit[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
    writeLines(nexus_edit,shuffled_alignment_path) # output the edited nexus file
    
    # Third, mask each alignment with the gaps and unknown characters from the original sequence 
    print("masking alignment")
    # for each alignment:
    #     - copy the sequence out from the new alignment: temp <- as.numeric(new_aln$X)
    #     - replace the non 18s with the generated sequence of the right length: temp[which(new_aln$X !=18)] <- new_seq
    #     - replace the new seq into the new aln: new_aln$X <- temp
    # Open the new alignment as a nexus file
    n_new <- read.nexus.data(shuffled_alignment_path)
    # Get the names of all the sequences
    seq_names <- names(n_new)
    print(paste0("number of names: ",length(seq_names)))
    print(paste0("number of unique names: ",length(unique(seq_names))))
    # Iterate through the names
    for (seq_name in seq_names){
      original_seq <- n[[seq_name]] # get the original empirical sequence
      new_seq <- n_new[[seq_name]] # get the new simulated sequence that has the same name
      gap_inds <- which(original_seq == "-") # find out which sites are a gap in the original alignment
      unknown_inds <- which(original_seq == "?") # find out which sites are unknown in the original alignment
      new_seq[gap_inds] <- "-" # add the gaps into the simulated alignment
      new_seq[unknown_inds] <- "?" # add the unknowns into the simulated alignment
      n_new[[seq_name]] <- new_seq
    }
    # Output the final alignment (same parameters and gaps as input alignment) as a nexus file
    print("output nexus file")
    write.nexus.data(n_new,file = bootstrap_alignment_path, format = "dna", interleaved = FALSE)
    # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
    nexus_edit <- readLines(bootstrap_alignment_path) # open the new nexus file
    ind <- grep("BEGIN DATA",nexus_edit)+2 # find which line
    nexus_edit[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
    writeLines(nexus_edit,bootstrap_alignment_path) # output the edited nexus file
    
  }
  # Run all the test statistics
  # bootstrap_id will be "bootstrapReplicateXXXX" where XXXX is a number
  print("run test statistics")
  empirical.runTS(alignment_path = bootstrap_alignment_path, program_paths = program_paths, bootstrap_id = bootstrap_id)
}



empirical.bootstraps.wrapper <- function(empirical_alignment_path, program_paths, number_of_replicates){
  print("in empirical.bootstraps.wrapper")
  print(empirical_alignment_path)
  # Create output file names, the name of the loci and the file path of the loci location
  collated_ts_file <- paste0(dirname(empirical_alignment_path),"/",gsub(".nex","",basename(empirical_alignment_path)),"_collatedBSReplicates.csv")
  p_value_file  <- paste0(dirname(empirical_alignment_path),"/",gsub(".nex","",basename(empirical_alignment_path)),"_pValues.csv")
  loci_name <- gsub(".nex","",basename(empirical_alignment_path))
  alignment_folder <- dirname(empirical_alignment_path)
  
  # Check that the original alignment ran ok
  if (file.exists(paste0(empirical_alignment_path,".iqtree")) == FALSE || file.exists(paste0(empirical_alignment_path,".treefile")) == FALSE || file.exists(paste0(empirical_alignment_path,".lmap.eps")) == FALSE){
    print("need to rerun IQ-Tree")
    n <- read.nexus.data(empirical_alignment_path)
    n_taxa <- length(n)
    # Run IQ-tree on the alignment (if it hasn't already been run), and get the likelihood mapping results
    print("run IQTree")
    call.IQTREE.quartet(program_paths[["IQTree"]],empirical_alignment_path,n_taxa)
  }
  
  # Calculate the test statistics if it hasn't already been done
  ts_file <- paste0(dirname(empirical_alignment_path),"/",gsub(".nex","",basename(empirical_alignment_path)),"_testStatistics.csv")
  if (file.exists(ts_file) == FALSE){
    print("run test statistics")
    empirical.runTS(empirical_alignment_path, program_paths, bootstrap_id = "alignment")
  }
  
  #Check that the test statistic file ran ok 
  if (file.exists(ts_file) == FALSE){
    print("need to rerun test statistics")
    empirical.runTS(empirical_alignment_path, program_paths, bootstrap_id = "alignment")
  }
  
  #Extract the parameters from the .iqtree log file.
  print("get simulation params")
  params <- get.simulation.parameters(paste0(empirical_alignment_path,".iqtree"))
  
  # Create the bootstrap ids (pad out to 4 digits) - should be "bootstrapReplicateXXXX" where XXXX is a number
  bootstrap_ids <- paste0("bootstrapReplicate",sprintf("%04d",1:number_of_replicates))
  
  # Run all the bootstrap ids that HAVEN'T already been run (e.g. in previous attempts) using lapply (feed info into do1.empirical.parametric.bootstrap)
  # If the alignment doesn't exist OR the test statistic csv doesnt exist, this indicates a complete run has not previously been done 
  # These bootstrap replicates will thus be calculates
  # This should save A BUNCH of time because it means if the test statistic file exists, you don't have to run Splitstree four times
  print("run bootstrap replicates")
  bs_als <- paste0(alignment_folder,"/",loci_name,"_",bootstrap_ids,"/",loci_name,"_",bootstrap_ids,".nex")
  ts_csvs <- paste0(alignment_folder,"/",loci_name,"_",bootstrap_ids,"/",loci_name,"_",bootstrap_ids,"_testStatistics.csv")
  missing_als <- bs_als[!file.exists(bs_als)]
  missing_testStatistics <- bs_als[!file.exists(ts_csvs)]
  # Collate the missing files and identify the alignments to rerun
  all_to_run <- unique(c(missing_als,missing_testStatistics))
  ids_to_run <- bootstrap_ids[which((bs_als %in% all_to_run))]
  print(paste0("Number of alignments to run = ",length(ids_to_run)))
  print("run all previously-unrun bootstraps")
  if(length(ids_to_run)>0){
    mclapply(ids_to_run, do1.empirical.parametric.bootstrap, empirical_alignment_path = empirical_alignment_path, alignment_params = params, program_paths = program_paths, mc.cores = 25)
  }

  # Before you can collate all the bootstrap files, you need to check every bootstrap ran and rerun the failed ones
  # Generate the names of each alignment, the test statistics csvs, the .iqtree files, the treefiles, the likelihood mapping files
  # Check which of these files are missing
  print("check for missing alignments")
  bs_als <- paste0(alignment_folder,"/",loci_name,"_",bootstrap_ids,"/",loci_name,"_",bootstrap_ids,".nex")
  ts_csvs <- paste0(alignment_folder,"/",loci_name,"_",bootstrap_ids,"/",loci_name,"_",bootstrap_ids,"_testStatistics.csv")
  missing_als <- bs_als[!file.exists(bs_als)]
  missing_iqtree <- bs_als[!file.exists(paste0(bs_als,".iqtree"))]
  missing_tree <- bs_als[!file.exists(paste0(bs_als,".treefile"))]
  missing_lmap <- bs_als[!file.exists(paste0(bs_als,".lmap.eps"))]
  missing_testStatistics <- bs_als[!file.exists(ts_csvs)]
  # Collate the missing files and identify the alignments to rerun
  all_missing <- unique(c(missing_als,missing_iqtree,missing_tree,missing_lmap,missing_testStatistics))
  als_to_rerun <- bootstrap_ids[which((bs_als %in% all_missing))]
  print(paste0("Number of missing alignments to rerun = ",length(als_to_rerun)))
  # Rerun the missing als
  if (length(als_to_rerun)>0){
    mclapply(als_to_rerun, do1.empirical.parametric.bootstrap, empirical_alignment_path = empirical_alignment_path, alignment_params = params, program_paths = program_paths, mc.cores = 25)
  }
  
  # collate the bootstrap info into 1 file
  print("collate bootstraps")
  p_value_df <- collate.bootstraps(directory = alignment_folder, file.name = "testStatistics", id = loci_name, output.file.name = collated_ts_file)
  # add the column with the bootstrap replicates and "alignment"
  new_bootstrap_ids <- p_value_df$bootstrap_replicate_id # copy col
  aln_id <- grep(new_bootstrap_ids[!grepl("bootstrapReplicate",new_bootstrap_ids)],new_bootstrap_ids) # get which element of col is the alignment
  new_bootstrap_ids[aln_id] <- "alignment"
  p_value_df$bootstrap_id <- new_bootstrap_ids
  
  # Calculate the p-values and add them to the original test statistic dataframe
  print("calculate p values")
  # Open the original test statistic file
  ts_df <- read.csv(ts_file)
  # Calculate the p_values of the variables of interest
  # Calculate the p-values for each test statistic
  ts_df$PHI_mean_sig <- calculate.p_value(p_value_df$PHI_mean, p_value_df$bootstrap_id)
  ts_df$PHI_observed_sig <- calculate.p_value(p_value_df$PHI_observed, p_value_df$bootstrap_id)
  ts_df$x3seq_sig <- calculate.p_value(p_value_df$X3SEQ_num_distinct_recombinant_sequences, p_value_df$bootstrap_id)
  ts_df$prop_resolved_quartets_sig <- calculate.p_value(p_value_df$prop_resolved_quartets, p_value_df$bootstrap_id)
  ts_df$splittable_percentage_sig <- calculate.p_value(p_value_df$splittable_percentage, p_value_df$bootstrap_id)
  ts_df$pdm_difference_sig <- calculate.p_value(p_value_df$pdm_difference, p_value_df$bootstrap_id)
  ts_df$pdm_average_sig <- calculate.p_value(p_value_df$pdm_average, p_value_df$bootstrap_id)
  ts_df$sd_untrimmed_sig <- calculate.p_value(p_value_df$split_decomposition_untrimmed, p_value_df$bootstrap_id)
  ts_df$nn_untrimmed_sig <- calculate.p_value(p_value_df$neighbour_net_untrimmed, p_value_df$bootstrap_id)
  ts_df$sd_trimmed_sig <- calculate.p_value(p_value_df$split_decomposition_trimmed, p_value_df$bootstrap_id)
  ts_df$nn_trimmed_sig <- calculate.p_value(p_value_df$neighbour_net_trimmed, p_value_df$bootstrap_id)
  ts_df$mean_delta_q_sig  <- calculate.p_value(p_value_df$mean_delta_q, p_value_df$bootstrap_id)
  ts_df$median_delta_q_sig <- calculate.p_value(p_value_df$median_delta_q, p_value_df$bootstrap_id)
  ts_df$mode_delta_q_sig <- calculate.p_value(p_value_df$mode_delta_q, p_value_df$bootstrap_id)
  # Output the p-values file
  write.csv(ts_df,file = p_value_file, row.names = FALSE)
}


fix.gammaCategory.siteNums <- function(df,num){
  print("fix gamma categories")
  # Quick function to randomly add/subtract to make sure that the number of sites in each gamma category is correct
  if (sum(df$cat_sites)==num){
    # If the sum of the num of sites = num of sites, return the df (all is good)
    return(df)
  }
  else {
    while (sum(df$cat_sites)!=num){
      # randomly sample a row
      row <- df[sample(nrow(df),1),]
      row_id <- which(df$category==row$category)
      if (sum(df$cat_sites)<num){
        row$cat_sites <- row$cat_sites + 1
      }
      if (sum(df$cat_sites)>num){
        row$cat_sites <- row$cat_sites - 1
      }
      # replace the row
      df[row_id,]<-row
    }
    # outside the while loop: therefor return the gamma categories matrix
    return(df)
  }
}
