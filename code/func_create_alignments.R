# Functions to create alignments for simulations with a coalescent or phylogenetic approach

library(TreeSim)
library(phytools)
library(seqinr)
library(ape)
library(phangorn)

# Create a function to simulate phylogenetic alignments using a number of parameters (birth rate, molecular rate, set trees)
# sites_from_tree_2 is the proportion of the SECOND tree that will be included (provide a single value between 0 and 1)
phylo.make1 <- function(output_folder, ntaxa, nsites, birth_rate = 0.5, tree_age = 1, mol_rate, mol_rate_sd = 0.1, sites_from_tree_2, id){
  # Randomly select a death rate using a uniform distribution from 0 to 99% of the birth rate
  death_rate = runif(1,min = 0, max = (0.99*birth_rate))
  # 1. Simulate a tree
  # simulate a birth-death tree on a fixed number of extant taxa
  # n = number of taxa, numbsim = # of trees to simulate, lambda = speciation rate [good default = 0.5 from Duchenne and Lanfear (2015)]
  # mu = extinction rate (default for this project = 0)
  tree_sim <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = birth_rate, mu = death_rate, frac = 1, age = tree_age, mrca = FALSE)[[1]] #simulate chronogram
  tree_sim$edge.length <- tree_sim$edge.length * (tree_age / max(branching.times(tree_sim))) # adjust branches to exactly equal tree age (adjusts for rounding errors)
  # 2. Simulate rate variation
  # Default for mol_rate and mol_rate_sd = 0.1 as in Duchenne and Lanfear (2015)
  phylo_sim <- tree_sim
  # simulate phylogram
  phylo_sim$edge.length <- tree_sim$edge.length * rlnorm(length(tree_sim$edge.length), meanlog = log(mol_rate), sdlog = mol_rate_sd) # adjust branch lengths - mol rate will control the tree depth in substitutions per site
  
  # Perform a single SPR move to get a new tree 
  phylo_sim_2 <- rSPR(phylo_sim, moves=1) # perform a single SPR move at random
  # to get SPR distance between the trees: SPR.dist()
  # spr_dist <- sprdist
  
  # Calculate how many sites of each tree will be needed 
  K <- sites_from_tree_2
  J <- 1 - K # proportion of first tree that will be included
  J_sites <- nsites * J # find number of sites to model from first tree
  K_sites <- nsites * K # find number of sites to model from second tree
  if ((J_sites+K_sites) < nsites){ 
    # if there are less sites then necessary, randomly add the missing sites to either tree one or tree 2
    add <- nsites - (J_sites+K_sites)
    rand <- sample(c("J","K"),1)
    if (rand == "K"){
      K_sites <- K_sites + add
    } else if (rand == "J"){
      J_sites <- J_sites + add
    }
  } else if ((J_sites+K_sites) > nsites){
    # if there are more sites then necessary, randomly subtract the missing sites from either tree one or tree 2
    subtract <- (J_sites+K_sites) - nsites
    rand <- sample(c("J","K"),1)
    if (rand == "K"){
      K_sites <- K_sites - subtract
    } else if (rand == "J"){
      J_sites <- J_sites - subtract
    }
  }
  # Simulate the DNA alignment
  if ((K_sites > 0) && (J_sites > 0)){
    # if sites are present from both trees, create an alignment for each tree and concatenate the alignments together
    dna_sim_1 <- simSeq(phylo_sim,l = J_sites) # simulate sites along the first tree
    dna_sim_2 <- simSeq(phylo_sim_2,l = K_sites) # simulate sites along the second tree (should be the smaller number - K is proportion of tree 2)
    dnabin_1 <- as.DNAbin(dna_sim_1) # convert to DNAbin
    dnabin_2 <- as.DNAbin(dna_sim_2) # convert to DNAbin
    dna_bin <- cbind(dnabin_1,dnabin_2, check.names = TRUE, fill.with.gaps = TRUE, quiet = FALSE) # concatenate the two alignments
    dna_sim <- as.phyDat(dna_bin) # convert the alignment to 
  } else if ((K_sites > 0) && (J_sites == 0)){
    dna_sim <- simSeq(phylo_sim_2,l = K_sites) # if no sites on tree1, simulate sites only along the second tree
  } else if ((K_sites == 0) && (J_sites > 0)) {
    dna_sim <- simSeq(phylo_sim,l = J_sites) # if no sites on tree2, simulate sites only along the first tree
  }
  # Output all the files
  # Make an output name for the nexus file
  output_name_template <- paste0(output_folder,"alignment.nexus") # create a name for the output file
  write.phyDat(dna_sim,file = output_name_template, format = "nexus",interleaved = TRUE, datablock = FALSE) # write the output as a nexus file
  # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
  nexus <- readLines(output_name_template) # open the new nexus file
  ind <- grep("BEGIN CHARACTERS",nexus)+2 # find which line
  nexus[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
  writeLines(nexus,output_name_template) # output the edited nexus file
  
  # Save the first tree and a picture of the first tree
  pdf(file = paste0(output_folder,"tree1.pdf"))
  plot.phylo(phylo_sim)
  dev.off()
  write.tree(phylo_sim, file = paste0(output_folder,"tree1.treefile"), tree.names = TRUE)
  # Save the second tree and a picture of the second tree
  pdf(file = paste0(output_folder,"tree2.pdf"))
  plot.phylo(phylo_sim_2)
  dev.off()
  write.tree(phylo_sim_2, file = paste0(output_folder,"tree2.treefile"), tree.names = TRUE)
  
  # output a text file with all the parameters
  output_name_template <- paste0(output_folder,"params.csv") # create a name for the output file 
  row <- c("phylogenetic",ntaxa,nsites,"NA","NA","NA",birth_rate,death_rate,tree_age,mol_rate,mol_rate_sd,J,K,id) # gather up all the variables
  names <- c("method","n_taxa","n_sites","internal_recombination","external_recombination","mutation_rate","birth_rate","death_rate","tree_age","mean_molecular_rate",
             "sd_molecular_rate","proportion_tree1","proportion_tree2","id") # gather up the var names
  df <- data.frame(matrix(nrow=0,ncol=14)) # make an empty dataframe
  df <- rbind(df,row) # attach the info to the empty df
  names(df) <- names # rename it so it's pretty and also actually helpful
  write.csv(df, file = output_name_template) # write the csv so you can use it later. 
}


# Wrapper function - feeds relevant row into the phylo.fixedtrees.run1sim function
phylo.fixedtrees.wrapper <- function(index, dataframe, program_paths, tree_folder){
  relevant_row <- dataframe[index,]
  phylo.fixedtrees.run1sim(relevant_row, program_paths, tree_folder)
}


# Function to run one entire simulation using a phylogenetic framework : create the alignment, run test statistics, and save 
phylo.fixedtrees.run1sim <- function(row, program_paths, tree_folder){
  # row needs to include: output folder, n_sites, tree_age, tree1, tree2, proportion_tree2, id, rep
  
  # Call the function to make the output folder name, alignment name, and results file name
  al_folder <- phylo.fixedtrees.output.folder(row)[1]
  al_file <- phylo.fixedtrees.output.folder(row)[2]
  results_file <- phylo.fixedtrees.output.folder(row)[3]
  
  # Check to see if the folder for the alignment exists and if it doesn't exist, create it
  if (dir.exists(al_folder)==FALSE){
    dir.create(al_folder)
  }
  
  
  # Extract values for creating the phylogenetic alignment from the input row (convert to numeric so can use the elements for ~ maths things ~)
  output_folder <- row[["output_folder"]]
  n_sites <- as.numeric(row[["n_sites"]])
  tree_age <- as.numeric(row[["tree_age"]])
  tree1_name <- row[["tree1"]]
  tree2_name <- row[["tree2"]]
  K <- as.numeric(row[["proportion_tree2"]])
  id <- paste0(row[["id"]],"_",row[["rep"]])
  
  # Open the trees 
  tree1 <- open.fixed.tree(tree1_name,tree_folder)
  tree2 <- open.fixed.tree(tree2_name, tree_folder)
  
  # Create the alignment (if it hasn't already been created)
  if (file.exists(paste0(al_folder,"alignment.nexus")) == FALSE) {
    phylo.fixedtrees.make1(al_folder, n_sites, tree_age, tree1, tree1_name, tree2, tree2_name, K, id)
  } 
  
  # The alignment now definitely exists. Now you can run IQ-tree on the alignment (if it hasn't already been run)
  n_taxa <- length(tree1$tip.label) # get the number of taxa
  if (file.exists(paste0(al_folder,"alignment.nexus.iqtree")) == FALSE){
    # Call the scf version of the IQtree function - this will result in scf being run
    num_scf_quartets <- choose(n_taxa,4)
    scf <- calculate.sCF(program_paths[["IQTree"]], al_file, n_taxa, num_threads = "1", num_scf_quartets)
  }
    
  # Set wd to alignment folder - means that 3seq and Phi files will be saved into the folder with their alignment
  setwd(al_folder)
  # Get paths to PhiPac, 3SEQ
  phi_path <- program_paths[["Phi"]] # get path to phipack executable
  seq_path <- program_paths[["3seq"]] # get path to 3seq executable
  # Note that Phi and 3Seq will only be run if they haven't already been run (checks for log files)
  filetype = tail(strsplit(al_file,"\\.")[[1]],n=1) # extract file format
  # run PHIPACK and 3seq (depending on the file format, will need to convert to fasta)
  if (filetype == "fasta"){
    # if the alignment is already in fasta format, run PhiPack through R
    if (file.exists(paste0(al_folder,"Phi.log")) == FALSE){
      phi_command <- paste0(phi_path," -f ",al_file, " -v") # assemble system command
      system(phi_command) #call phipack
    }
    
    if (file.exists(paste0(al_folder,"3s.log")) == FALSE){
      seq_command <- paste0(seq_path," -f ", al_file)
      system(seq_command) #call 3SEQ
    }
  } else if (filetype == "nexus"){
    # Phipack only reads in Phylip or fasta format - need to convert if the alignment is a nexus file
    data = read.nexus.data(al_file) # read in nexus format alignment
    fasta.name <- paste0(al_file,".fasta") # make a name for the fasta alignment by adding .fasta (super original ;) )
    write.fasta(sequences = data,names = names(data), file.out = fasta.name) # output alignment as a fasta format
    if (file.exists(paste0(al_folder,"Phi.log")) == FALSE){
      phi_command <- paste0(phi_path," -f ",fasta.name, " -v") # assemble system command as above
      system(phi_command) # run PHI test on the new fasta alignment
    }
    
    if (file.exists(paste0(al_folder,"3s.log")) == FALSE){
      seq_command <- paste0(seq_path," -f ", fasta.name)
      system(seq_command) #call 3SEQ
    }
  }
  # Extract significance from Phi Pack output
  phi_file <- paste0(al_folder,"Phi.log")
  phi_file <- readLines(phi_file)
  ind      <- grep("p-Value",phi_file)
  phi_sig <- as.numeric(strsplit(phi_file[ind+3],":")[[1]][2])
  ind      <- grep("PHI Values",phi_file)
  phi_mean <- as.numeric(strsplit(phi_file[ind+4],"      ","")[[1]][2])
  phi_var <- as.numeric(strsplit(phi_file[ind+5],"      ","")[[1]][2])
  phi_obs <- as.numeric(strsplit(phi_file[ind+6],"      ","")[[1]][2])
  
  # Extract results output from 3Seq output
  seq_file <- paste0(al_folder,"3s.log")
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
  
  # Extract quartet mapping
  iq_log_path <- paste0(al_file,".iqtree")
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
  pdmm <- as.matrix(mldist.pdm(al_file)) # open pairwise distance matrix as a matrix
  deltaplot_results <- delta.plot(pdmm, k = 51, plot = FALSE) # calculate the delta.plot
  counts <- deltaplot_results$counts
  intervals <- seq(0,1,(1/(length(counts)-1)))
  deltaplot_df <- data.frame(intervals,counts)
  names(deltaplot_df) <- c("intervals","counts")
  deltaplot_df_name <- phylo.fixedtrees.output.folder(row)[4]
  write.csv(deltaplot_df,file = deltaplot_df_name)
  # Want to calculate the mean and median delta q value - unfortunately the delta.plot function doesn't output raw data, so make a pseudo data set using the histogram values
  mean_dq <- mean(rep(deltaplot_df$intervals,deltaplot_df$counts)) # turn the interval data into a long list of "raw" values and calculate the mean
  median_dq <- median(rep(deltaplot_df$intervals,deltaplot_df$counts)) # turn the interval data into a long list of "raw" values and calculate the median
  mode_dq <- deltaplot_df[order(deltaplot_df$count, decreasing = TRUE),][1,1] # sort the dataframe by count values and extract the mode
  
  # Run my test statistics
  # run pdm ratio (TS1) (modified splittable percentage)
  splittable_percentage <- pdm.ratio(iqpath = program_paths[["IQTree"]], path = al_file)
  # run normalised.pdm.difference.sum (TS2a) (sum of difference of normalised matrix)
  npds <- normalised.pdm.diff.sum(iqpath = program_paths[["IQTree"]], path = al_file)
  # run normalised pdm difference average (TS2b) (mean of difference of normalised matrix)
  npdm <- normalised.pdm.diff.mean(iqpath = program_paths[["IQTree"]], path = al_file)
  # Run trimmed and untrimmed versions of the split decomposition and NeighborNet tree proportion
  sd_untrimmed <- tree.proportion(iqpath = program_paths[["IQTree"]], splitstree_path = program_paths[["SplitsTree"]], path = al_file, network_algorithm = "split decomposition", trimmed = FALSE)
  nn_untrimmed <- tree.proportion(iqpath = program_paths[["IQTree"]], splitstree_path = program_paths[["SplitsTree"]], path = al_file, network_algorithm = "neighbournet", trimmed = FALSE)
  sd_trimmed <- tree.proportion(iqpath = program_paths[["IQTree"]], splitstree_path = program_paths[["SplitsTree"]], path = al_file, network_algorithm = "split decomposition", trimmed = TRUE)
  nn_trimmed <- tree.proportion(iqpath = program_paths[["IQTree"]], splitstree_path = program_paths[["SplitsTree"]], path = al_file, network_algorithm = "neighbournet", trimmed = TRUE)
  
  # Collect results
  all_files <- list.files(al_folder) # get a list of all the files
  ind <- grep("params",all_files) # find which of those files is the parameters file
  params_csv <- read.csv(all_files[ind]) # open the parameter file
  proportion_tree_1 <- params_csv$proportion_tree1 # extract proportion of alignment from tree 1
  n_taxa <- length(tree1$tip.label) #ntaxa = number of labels on the tree
  # Make somewhere to store the results
  df_names <- c("alignment","method","n_taxa","n_sites","tree_age","tree1","proportion_tree1","tree2","proportion_tree2","id",
                "PHI_mean","PHI_variance","PHI_observed","PHI_sig","3SEQ_num_recombinant_triplets","3SEQ_num_distinct_recombinant_sequences","3SEQ_p_value","num_quartets",
                "num_resolved_quartets","prop_resolved_quartets","num_partially_resolved_quartets","num_unresolved_quartets", "splittable_percentage","pdm_difference",
                "pdm_average","split_decomposition_untrimmed", "neighbour_net_untrimmed", "split_decomposition_trimmed","neighbour_net_trimmed","mean_delta_q","median_delta_q",
                "mode_delta_q", "sCF_mean", "sCF_median")
  df <- data.frame(matrix(nrow=0,ncol=length(df_names))) # create an empty dataframe of the correct size
  # collect all the variables
  op_row <- c(al_file, "phylogenetic_fixedTrees" ,n_taxa,row[["n_sites"]], row[["tree_age"]], tree1_name,proportion_tree_1,
              tree2_name, row[["proportion_tree2"]], paste0(row[["id"]],"_",row[["rep"]]), phi_mean, phi_var, phi_obs, phi_sig, 
              num_trips, num_dis, seq_sig, total_q, resolved_q, prop_resolved, partly_resolved_q, unresolved_q, splittable_percentage, 
              npds, npdm, sd_untrimmed, nn_untrimmed, sd_trimmed, nn_trimmed, mean_dq, median_dq, mode_dq, sCF$mean_scf, sCF$median_scf)
  df <- rbind(df,op_row,stringsAsFactors = FALSE) # place row in dataframe
  names(df) <- df_names # add names to the df so you know what's what
  write.csv(df,file = results_file)
}


# Name output folders and files
phylo.fixedtrees.output.folder <- function(row){
  # Extract values for creating the phylogenetic alignment from the input row (convert to numeric so can use the elements for ~ maths things ~)
  n_sites <- as.numeric(row[["n_sites"]])
  tree_age <- as.numeric(row[["tree_age"]])
  tree1_str <- gsub("_","-",row[["tree1"]])
  tree2_str <- gsub("_","-",row[["tree2"]])
  K <- as.numeric(row[["proportion_tree2"]])
  id <- paste0(row[["id"]],"_",row[["rep"]])
  # Create an output folder name using the variables
  output_folder <- paste0(row[["output_folder"]],"Phylo_FixedTrees_",n_sites,"_",tree_age,"_",tree1_str,"_",tree2_str,"_",K,"_",id,"/")
  nexus_file <- paste0(output_folder,"alignment.nexus")
  output_file <- paste0(output_folder,"testStatistics.csv")
  deltaplot_file <- paste0(output_folder,"deltaplot_histogram.csv")
  return(c(output_folder,nexus_file,output_file,deltaplot_file))
}

open.fixed.tree <- function(tree_name,tree_folder){
  tree_path <- paste0(tree_folder, tree_name, ".txt")
  phylogenetic_tree <- read.tree(tree_path)
  return(phylogenetic_tree)
}

# Create a function to make phylogenetic alignments (as outlined in simulation scheme)
# K is the proportion of the SECOND tree that will be included (provide a single value)
phylo.fixedtrees.make1 <- function(output_folder, n_sites, tree_age, tree1, tree1_str, tree2, tree2_str, K, id){
  # Scale tree1 and tree2 to have total length of tree_age
  tree1$edge.length <- tree1$edge.length * (tree_age / max(branching.times(tree1))) # adjust branches to have maximum branching length exactly equal to tree age
  tree2$edge.length <- tree2$edge.length * (tree_age / max(branching.times(tree2))) # adjust branches to have maximum branching length exactly equal to tree age

  # Calculate how many sites of each tree will be needed 
  J <- 1 - K # proportion of first tree that will be included
  J_sites <- n_sites * J # find number of sites to model from first tree
  K_sites <- n_sites * K # find number of sites to model from second tree
  if ((J_sites+K_sites) < n_sites){ 
    # if there are less sites then necessary, randomly add the missing sites to either tree one or tree 2
    add <- n_sites - (J_sites+K_sites)
    rand <- sample(c("J","K"),1)
    if (rand == "K"){
      K_sites <- K_sites + add
    } else if (rand == "J"){
      J_sites <- J_sites + add
    }
  } else if ((J_sites+K_sites) > n_sites){
    # if there are more sites then necessary, randomly subtract the missing sites from either tree one or tree 2
    subtract <- (J_sites+K_sites) - n_sites
    rand <- sample(c("J","K"),1)
    if (rand == "K"){
      K_sites <- K_sites - subtract
    } else if (rand == "J"){
      J_sites <- J_sites - subtract
    }
  }
  # Simulate the DNA alignment
  if ((K_sites > 0) && (J_sites > 0)){
    # if sites are present from both trees, create an alignment for each tree and concatenate the alignments together
    dna_sim_1 <- simSeq(tree1,l = J_sites) # simulate sites along the first tree
    dna_sim_2 <- simSeq(tree2,l = K_sites) # simulate sites along the second tree (should be the smaller number - K is proportion of tree 2)
    dnabin_1 <- as.DNAbin(dna_sim_1) # convert to DNAbin
    dnabin_2 <- as.DNAbin(dna_sim_2) # convert to DNAbin
    dna_bin <- cbind(dnabin_1,dnabin_2, check.names = TRUE, fill.with.gaps = TRUE, quiet = FALSE) # concatenate the two alignments
    dna_sim <- as.phyDat(dna_bin) # convert the alignment to phydat 
  } else if ((K_sites > 0) && (J_sites == 0)){
    dna_sim <- simSeq(tree2,l = K_sites) # if no sites on tree1, simulate sites only along the second tree
  } else if ((K_sites == 0) && (J_sites > 0)) {
    dna_sim <- simSeq(tree1,l = J_sites) # if no sites on tree2, simulate sites only along the first tree
  }
  # Output all the files
  # Make an output name for the nexus file
  output_name_template <- paste0(output_folder,"alignment.nexus") # create a name for the output file
  write.phyDat(dna_sim,file = output_name_template, format = "nexus",interleaved = TRUE, datablock = FALSE) # write the output as a nexus file
  # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
  nexus <- readLines(output_name_template) # open the new nexus file
  ind <- grep("BEGIN CHARACTERS",nexus)+2 # find which line
  nexus[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
  writeLines(nexus,output_name_template) # output the edited nexus file
  
  # Save the first tree and a picture of the first tree
  pdf(file = paste0(output_folder,"tree1.pdf"))
  plot.phylo(tree1)
  dev.off()
  write.tree(tree1, file = paste0(output_folder,"tree1.treefile"), tree.names = TRUE)
  # Save the second tree and a picture of the second tree
  pdf(file = paste0(output_folder,"tree2.pdf"))
  plot.phylo(tree2)
  dev.off()
  write.tree(tree2, file = paste0(output_folder,"tree2.treefile"), tree.names = TRUE)
  
  # Make the parameters file
  n_taxa <- length(tree1$tip.label) # get the number of taxa for the parameters csv
  # output a text file with all the parameters
  output_name_template <- paste0(output_folder,"params.csv") # create a name for the output file 
  row <- c("phylogenetic_fixedTrees",n_taxa,n_sites,tree_age,tree1_str,J,tree2_str,K,id) # gather up all the variables
  names <- c("method","n_taxa","n_sites","tree_age","tree1","proportion_tree1","tree2","proportion_tree2","id") # gather up the var names
  df <- data.frame(matrix(nrow=0,ncol=9)) # make an empty dataframe
  df <- rbind(df,row) # attach the info to the empty df
  names(df) <- names # rename it so it's pretty and also actually helpful
  write.csv(df, file = output_name_template) # write the csv so you can use it later. 
}



