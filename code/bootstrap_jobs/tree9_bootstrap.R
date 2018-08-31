# R code to create dataframes with a row for each simulation to run, simulate those alignments and output a csv with results of all test statistics
# Bootstrap - run tree 9

# Open packages
library(parallel)
library(seqinr)
library(ape)
library(TreeSim)
library(phytools)
library(phangorn)
library(base)
library(ggplot2)
library(reshape2)
# library(tictoc) # library for measuring timings! tic("label") to start, toc() to stop

# Remember to have downloaded and tested the following programs: SplitsTree4, SimBac, IQ-Tree, PhiPack and 3Seq. 
# 3Seq must have been associated with a P-value table for it to run properly
# If it's not working on soma, you (Caitlin) might have forgotten about the LD_LIBRARY_PATH again: run these two lines of code
# LD_LIBRARY_PATH=$LD_LIBRARY_PATH:usr/local/lib:usr/lib/x86_64-linux-gnu
# export LD_LIBRARY_PATH

# Set parameters
run_location <- "soma"
op_folder <- "/data/caitlin/treelikeness/output/"
results_folder <- "/data/caitlin/treelikeness/results/"
maindir <- "/data/caitlin/treelikeness/"
exec_paths <- c("/data/caitlin/linux_executables/3seq/3seq","/data/caitlin/linux_executables/iqtree/bin/iqtree","/data/caitlin/linux_executables/PhiPack/Phi",
                "/data/caitlin/linux_executables/SimBac/SimBac","/data/caitlin/splitstree4/SplitsTree")
names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")
network_functions <- "code/func_split_decomposition.R"
run_id <- "soma_fullSet"

# Set working directory
setwd(maindir)
# Source files for functions
source(paste0(maindir,network_functions))
source(paste0(maindir,"code/func_test_statistic.R"))
source(paste0(maindir,"code/func_create_alignments.R"))
source(paste0(maindir,"code/func_process_data.R"))
source(paste0(maindir,"code/func_parametric_bootstrap.R"))
tree_folder <- paste0(maindir,"trees/")

## For fourth set of plots:
# Make empty dataframe:
plot4_df <- data.frame((matrix(ncol = 8, nrow = 0)))
names(plot4_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
# Parameters that are the same for each set of trees:
output_folder <- op_folder
n_sites <- 1300
tree_age <- c(0.05, 0.1, 0.5, 1)
proportion_tree2 <- seq(0,0.5,0.1)
id <- "plot4tree9"
rep <- 1:100
# tree1_vector <- c("08taxa_balanced_LHS","08taxa_balanced_LHS",
#                   "08taxa_intermediate_LHS","08taxa_intermediate_LHS",
#                   "08taxa_unbalanced_LHS","08taxa_unbalanced_LHS",
#                   "08taxa_balanced_LHS","08taxa_intermediate_LHS","08taxa_unbalanced_LHS")
# tree2_vector <- c("08taxa_balanced_RHS_reciprocal_close_1event","08taxa_balanced_RHS_nonreciprocal_close_1event",
#                   "08taxa_intermediate_RHS_reciprocal_close_1event","08taxa_intermediate_RHS_nonreciprocal_close_1event",
#                   "08taxa_unbalanced_RHS_reciprocal_close_1event","08taxa_unbalanced_RHS_nonreciprocal_close_1event",
#                   "08taxa_balanced_LHS","08taxa_intermediate_LHS","08taxa_unbalanced_LHS")

tree1_vector <- c("08taxa_unbalanced_LHS")
tree2_vector <- c("08taxa_unbalanced_LHS")

tree_id <- 1:length(tree1_vector)
for (i in tree_id){
  tree1_temp <- tree1_vector[[i]]
  tree2_temp <- tree2_vector[[i]]
  temp_df <- expand.grid(op_folder,n_sites,tree_age,tree1_temp,tree2_temp,proportion_tree2,id,rep, stringsAsFactors = FALSE)
  names(temp_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
  plot4_df <- rbind(plot4_df,temp_df, stringsAsFactors = FALSE)
}
mclapply(1:nrow(plot4_df), phylo.fixedtrees.wrapper, plot4_df, exec_paths, tree_folder, mc.cores = 4) # mclapply for phylo with fixed trees
# Collect the folders that contain the alignments for plot4
all_folders <- list.dirs(op_folder, recursive = FALSE, full.names = TRUE) # get all the directory names in the output folder
inds <- grep(id,all_folders) # find which indexes the plot4 (bootstrap) folders are at 
plot4_folders <- all_folders[inds] # get the bootstrap folders
plot4_folders <- paste0(plot4_folders,"/") # add the slash to the end so it's a path to the directory. Bootstrap function just adds "alignment.nex" not the slash.
plot4_toRun <- c() # create an empty list to store the folders that need the bootstrap run
# Check each simulation from the plot4_df for a p value csv
for (folder in plot4_folders){
  p_file <- paste0(folder,"p_value.csv")
  if (file.exists(p_file) == FALSE) {
    # if there's no p value csv, there's no bootstrap: add to the list of bootstraos to run
    plot4_toRun <- c(plot4_toRun, folder)
  }
}
# Apply the parametric bootstrap function to the folders without a bootstrap
mclapply(plot4_toRun, phylo.parametric.bootstrap, 199, exec_paths[["IQTree"]], exec_paths[["SplitsTree"]], exec_paths[["Phi"]], exec_paths[["3seq"]], mc.cores = 4) # run all the bootstraps!

