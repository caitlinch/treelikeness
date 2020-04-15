# R code to create dataframes with a row for each simulation to run, simulate those alignments and output a csv with results of all test statistics
# Sourcing this file will run all of the simulations below
# Final result is a seires of folder, one per simulation, each containing a simulated alignment, an IQ-Tree file tree and log, a csv containing all
# parameters from the simulation, a csv containing test statistics calculated for that alignment and a number of files generated from the programs used.
# For the third experiment, a number of parametric bootstrap replicates will also be generated, each nested inside the alignment folder in a separate folder.



# Remember to have downloaded and tested the following programs: SplitsTree4, IQ-Tree, PhiPack and 3Seq. 
# 3Seq must have been associated with a P-value table for it to run properly
# IQ-Tree version must be 1.7 or above to run site concordance factors



##### Step 1: Open packages #####
library(parallel)
library(seqinr)
library(ape)
library(TreeSim)
library(phytools)
library(phangorn)



##### Step 2: Uncomment and set the file paths for output folders, executables, and the identifying name for this run #####
# op_folder <- the folder where simulated alignments and output from analysis (e.g. IQ-Tree output files, 3seq output files) will be placed
# results_folder <- the folder where the result csvs will be placed. I use same results folder for parts 1-4.
# maindir <- "treelikeness" repository location
# run_id <- the key for this simulation. Will be in the names for outputs and the results csvs. Must NOT contain any underscores ("_")
# num_cores <- the number of cores to use. 1 for a single core (wholly sequential), or higher if using parallelisation.
# num_reps <- the number of of bootstrap replicates to perform (e.g if num_reps = 199, for each alignment 199 new alignments will be
#             simulated with the same parameters and used to determine whether the treelikeness of the original alignment is expected or not)
# max_reps <- max_reps is the number of simulated alignments for each set of simulation parameters that will have the parametric bootstrap applied
#             e.g. if max_reps is 20, 20 alignments for each set of simulation parameters will be fed into the parametric bootstrap calculation function
# exec_paths <- location to each executable within the folder. Attach the names of the executables so the paths can be accessed by name. 
#               DO NOT change the order of the executables.
# The SplitsTree executable path can be tricky to find: 
#       - in MacOS, the path is "SplitsTree.app/Contents/MacOS/JavaApplicationStub" (assuming you are in the same directory as the application)
#       - in Linux, after installing and navigating into the folder it's simply "SplitsTree"

# op_folder <- "" 
# results_folder <- ""
# maindir <- ""
# run_id <- ""
# num_cores <- 30
# num_reps <- 199
# max_reps <- 20

# Create a vector with all of the executable file paths
# The executables must be in the following order: 3seq, IQ-Tree, PHIPack, SplitsTree
# if executables are in different folders, remove the "exec_folder <- '' " line and the "exec_paths <- paste()" line and provide an individual path to each executable within exec_paths
# To access a path: exec_paths[["name"]]

# exec_folder <- "/path/to/executables/folder/"
# exec_paths <- c("3seq_executable","IQ-Tree_executable","PHIpack_executable","SplitsTree_executable") # MUST be in this order
# exec_paths <- paste0(exec_folder,exec_paths)


#__________________________________________Caitlin's paths (delete these if you're not Caitlin)______________________________________
op_folder <- "/data/caitlin/treelikeness/output_20200304/"
results_folder <- "/data/caitlin/treelikeness/results_20200304/"
maindir <- "/data/caitlin/treelikeness/"
run_id <- "sCF"
num_cores <- 20
num_reps <- 199
max_reps <- 20

# Create a vector with all of the executable file paths
# To access a path: exec_paths[["name"]]
exec_paths <- c("/data/caitlin/linux_executables/3seq/3seq","/data/caitlin/linux_executables/iqtree-2.0-rc1-Linux/bin/iqtree",
                "/data/caitlin/linux_executables/PhiPack/Phi","/data/caitlin/splitstree4/SplitsTree")
#____________________________________________________________________________________________________________________________________



##### Step 3: Source function files #####
# Set working directory
setwd(maindir)

# Attach the names used throughout the code and functions to the exec_paths object
# Do not edit these names - it will cause internal functions and thus calculations to fail
names(exec_paths) <- c("3seq","IQTree","Phi","SplitsTree")

# Source files for functions
source(paste0(maindir,"code/func_test_statistic.R"))
source(paste0(maindir,"code/func_create_alignments.R"))
source(paste0(maindir,"code/func_process_data.R"))
source(paste0(maindir,"code/func_parametric_bootstrap.R"))
tree_folder <- paste0(maindir,"trees/")




##### Step 4: Create parameter dataframes for the set of simulations and run them #####
### Create dataframe for the final set of simulations (fixed trees)
### Each row needs to include: output_folder, n_sites, tree_age, mean_molecular_rate, sd_molecular_rate, tree1, tree2, proportion_tree2,id,rep

## For first experiment:
# - What effect does varying the type of introgression event have on the treelikeness score?

# Make empty dataframe:
exp1_df <- data.frame((matrix(ncol = 8, nrow = 0)))
names(exp1_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
# Parameters that are the same for each set of trees:
output_folder <- op_folder
n_sites <- 1300
tree_age <- c(0.05, 0.1, 0.5, 1)
proportion_tree2 <- 0.5
id <- "exp1"
rep <- 1:100
tree1_vector <- c("08taxa_balanced_LHS","08taxa_balanced_LHS","08taxa_balanced_LHS","08taxa_balanced_LHS","08taxa_balanced_LHS","08taxa_balanced_LHS",
                  "08taxa_balanced_LHS")
tree2_vector <- c("08taxa_balanced_RHS_reciprocal_close_1event","08taxa_balanced_RHS_reciprocal_divergent_1event","08taxa_balanced_RHS_reciprocal_ancient_1event",
                  "08taxa_balanced_RHS_nonreciprocal_close_1event","08taxa_balanced_RHS_nonreciprocal_divergent_1event","08taxa_balanced_RHS_nonreciprocal_ancient_1event",
                  "08taxa_balanced_LHS")
# Expand parameters into dataframe
tree_id <- 1:7
for (i in tree_id){
  tree1_temp <- tree1_vector[[i]]
  tree2_temp <- tree2_vector[[i]]
  temp_df <- expand.grid(op_folder,n_sites,tree_age,tree1_temp,tree2_temp,proportion_tree2,id,rep, stringsAsFactors = FALSE)
  names(temp_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
  exp1_df <- rbind(exp1_df,temp_df, stringsAsFactors = FALSE)
}
# Run simulation
mclapply(1:nrow(exp1_df), phylo.fixedtrees.wrapper, exp1_df, exec_paths, tree_folder, mc.cores = num_cores) # mclapply for phylo with fixed trees


## For second experiment:
# - What effect on the treelikeness score does increasing the number of introgression events have?
#       - Fix the number of sites
#       - Starting with a balanced 8-taxon tree, perform 0 - 8 simultaneous close introgression events
#         (for both reciprocal and non-reciprocal events)
#       - Fix the proportion of tree 2 at 50% and the proportion of tree 1 at 50%
#       - Vary the tree depth from 0.05 to 1 substitution per site
#       - Perform 100 replicates for each set of parameters

# Make empty dataframe:
exp2_df <- data.frame((matrix(ncol = 8, nrow = 0)))
names(exp2_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
# Parameters that are the same for each set of trees:
output_folder <- op_folder
n_sites <- 1300
tree_age <- c(0.05, 0.1, 0.5, 1)
proportion_tree2 <- (0.5)
id <- "exp2"
rep <- 1:100
tree1_vector <- c("32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS",
                  "32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS",
                  "32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS",
                  "32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS",
                  "32taxa_balanced_LHS")
tree2_vector <- c("32taxa_balanced_RHS_reciprocal_close_1event","32taxa_balanced_RHS_reciprocal_close_2event",
                  "32taxa_balanced_RHS_reciprocal_close_3event","32taxa_balanced_RHS_reciprocal_close_4event",
                  "32taxa_balanced_RHS_reciprocal_close_5event","32taxa_balanced_RHS_reciprocal_close_6event",
                  "32taxa_balanced_RHS_reciprocal_close_7event","32taxa_balanced_RHS_reciprocal_close_8event",
                  "32taxa_balanced_RHS_nonreciprocal_close_1event","32taxa_balanced_RHS_nonreciprocal_close_2event",
                  "32taxa_balanced_RHS_nonreciprocal_close_3event","32taxa_balanced_RHS_nonreciprocal_close_4event",
                  "32taxa_balanced_RHS_nonreciprocal_close_5event","32taxa_balanced_RHS_nonreciprocal_close_6event",
                  "32taxa_balanced_RHS_nonreciprocal_close_7event","32taxa_balanced_RHS_nonreciprocal_close_8event",
                  "32taxa_balanced_LHS")
# Expand parameters into dataframe
tree_id <- 1:17
for (i in tree_id){
  tree1_temp <- tree1_vector[[i]]
  tree2_temp <- tree2_vector[[i]]
  temp_df <- expand.grid(op_folder,n_sites,tree_age,tree1_temp,tree2_temp,proportion_tree2,id,rep, stringsAsFactors = FALSE)
  names(temp_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
  exp2_df <- rbind(exp2_df,temp_df, stringsAsFactors = FALSE)
}
# Run simulation
mclapply(1:nrow(exp2_df), phylo.fixedtrees.wrapper, exp2_df, exec_paths, tree_folder, mc.cores = num_cores) # mclapply for phylo with fixed trees

# Collect the folders that contain the relevant alignments for second experiment and run the parametric bootstraps.
# Can't run all 100 alignments with 199 parametric bootstraps so subset this and only run max_reps replicates per parameter combination
# Need to filter all the alignments created for the analysis to just get the ones used for the parametric bootstrap - rep = 1:max_reps inclusive
all_folders <- paste0(op_folder,list.dirs(op_folder, recursive = FALSE, full.names = FALSE)) # get all the directory names in the output folder
inds <- grep(id,all_folders) # find which indexes the exp3 (bootstrap) folders are at 
exp2_folders <- all_folders[inds] # get the bootstrap folders
# Extract the exp2 alignments with a replicate number of max_reps or less
exp2_folders_bs <- unlist(lapply(exp2_folders, reject.high.reps, max_reps))
# Format the selected alignments so the parametric bootstrap can be carried out
exp2_folders_bs <- paste0(exp2_folders_bs,"/") # add the slash to the end so it's a path to the directory. Bootstrap function just adds "alignment.nex" not the slash.
exp2_toRun <- c() # create an empty list to store the folders that need the bootstrap run
# Check each simulation from the exp3_df for a p value csv
for (folder in exp2_folders_bs){
  p_file <- paste0(folder,"p_value.csv")
  if (file.exists(p_file) == FALSE) {
    # if there's no p value csv, there's no bootstrap: add to the list of bootstraps to run
    exp2_toRun <- c(exp2_toRun, folder)
  }
}
# Apply the parametric bootstrap function
mclapply(exp2_toRun, phylo.parametric.bootstrap, n_reps = num_reps, exec_paths[["IQTree"]], exec_paths[["SplitsTree"]], exec_paths[["Phi"]], exec_paths[["3seq"]], mc.cores = num_cores) # run all the bootstraps!


## For third experiment:

# - What effect does varying tree depth have on the treelikeness score? What effect does increasing the proportion of introgressed DNA have
#   on the treelikeness score?
#       - Fix the number of sites
#       - Starting with a balanced 8-taxon tree, perform introgression event (events may be reciprocal or non-reciprocal.
#         Events vary in location in the tree and may be close, divergent or ancient introgression.)
#       - Vary the proportion of tree 2 from 0 - 50% in 1% increments (meaning initial proportion of tree 1 is 100% and final
#         proportion of tree 1 is 50%)
#       - Vary the tree depth from 0.05 to 1 substitution per site
#       - Perform 100 replicates for each set of parameters
# - Investigating performance of test statistics using a parametric bootstrap
#       - Fix the number of sites
#       - Fix the tree as a balanced 8 taxon tree with one close introgression event (either reciprocal or non-reciprocal).
#       - Vary the proportion of tree 2 from 0 - 50% in 1% increments (meaning initial proportion of tree 1 is 100% and final
#         proportion of tree 1 is 50%)
#       - Vary the tree depth from 0.05 to 1 substitution per site
#       - Perform 100 replicates for each set of parameters
#       - Perform 199 parametric bootstraps for each alignment where proportion of tree 2 = 0%, 10%, 20%, 30%, 40% and 50%

# Make empty dataframe:
exp3_df <- data.frame((matrix(ncol = 8, nrow = 0)))
names(exp3_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
output_folder <- op_folder
n_sites <- 1300
tree1_vector <- c("08taxa_balanced_LHS","08taxa_balanced_LHS","08taxa_balanced_LHS")
tree2_vector <- c("08taxa_balanced_RHS_reciprocal_close_1event","08taxa_balanced_RHS_nonreciprocal_close_1event","08taxa_balanced_LHS")
tree_age <- c(0.05, 0.1, 0.5, 1)
proportion_tree2 <- seq(0,0.5,0.01)
id <- "exp3"
rep <- 1:100
# Expand parameters into dataframe
tree_id <- 1:3
for (i in tree_id){
  tree1_temp <- tree1_vector[[i]]
  tree2_temp <- tree2_vector[[i]]
  temp_df <- expand.grid(op_folder,n_sites,tree_age,tree1_temp,tree2_temp,proportion_tree2,id,rep, stringsAsFactors = FALSE)
  names(temp_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
  exp3_df <- rbind(exp3_df,temp_df, stringsAsFactors = FALSE)
}
# Run simulation
mclapply(1:nrow(exp3_df), phylo.fixedtrees.wrapper, exp3_df, exec_paths, tree_folder, mc.cores = num_cores) # mclapply for phylo with fixed trees

# Collect the folders that contain the alignments for third experiment and run the parametric bootstraps. 
all_folders <- paste0(op_folder,list.dirs(op_folder, recursive = FALSE, full.names = FALSE)) # get all the directory names in the output folder
inds <- grep(id,all_folders) # find which indexes the exp3 (bootstrap) folders are at 
exp3_folders <- all_folders[inds] # get the bootstrap folders
exp3_folders <- paste0(exp3_folders,"/") # add the slash to the end so it's a path to the directory. Bootstrap function just adds "alignment.nex" not the slash.
# extract all folders you want a bootstrap for (i.e. proportion_tree2 in seq(0,0.5,0.1) e.g. 10%)
trimmed_exp3_folders <- unlist(lapply(exp3_folders, reject.exp3.bootstrap.run))
trimmed_exp3_folders <- unlist(lapply(trimmed_exp3_folders, reject.high.reps, max_reps))
exp3_toRun <- c() # create an empty list to store the folders that need the bootstrap run
# Check each simulation from the exp3_df for a p value csv
for (folder in trimmed_exp3_folders){
  p_file <- paste0(folder,"p_value.csv")
  if (file.exists(p_file) == FALSE) {
    # if there's no p value csv, there's no bootstrap: add to the list of bootstraps to run
    exp3_toRun <- c(exp3_toRun, folder)
  }
}
# Apply the parametric bootstrap function
mclapply(exp3_toRun, phylo.parametric.bootstrap, n_reps = num_reps, exec_paths[["IQTree"]], exec_paths[["SplitsTree"]], exec_paths[["Phi"]], exec_paths[["3seq"]], mc.cores = num_cores) # run all the bootstraps!



##### Step 5: Save the parameter dataframes #####
op_name <- paste0(results_folder,"exp1_input_parameters_",run_id,".csv")
write.csv(exp1_df,file=op_name)
op_name <- paste0(results_folder,"exp2_input_parameters_",run_id,".csv")
write.csv(exp2_df,file=op_name)
op_name <- paste0(results_folder,"exp3_input_parameters_",run_id,".csv")
write.csv(exp3_df,file=op_name)



