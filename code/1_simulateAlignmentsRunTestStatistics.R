# R code to create dataframes with a row for each simulation to run, simulate those alignments and output a csv with results of all test statistics

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
library(tictoc) # library for measuring timings! tic("label") to start, toc() to stop

# Remember to have downloaded and tested the following programs: SplitsTree4, SimBac, IQ-Tree, PhiPack and 3Seq. 
# 3Seq must have been associated with a P-value table for it to run properly
# If it's not working on soma, you (Caitlin) might have forgotten about the LD_LIBRARY_PATH again: run these two lines of code
# LD_LIBRARY_PATH=$LD_LIBRARY_PATH:usr/local/lib:usr/lib/x86_64-linux-gnu
# export LD_LIBRARY_PATH

# run_location <- "mac"
# run_location <- "nci"
run_location <- "soma"

times <- c()
time_ids <- c()

tic("initialising")
if (run_location == "mac"){
  	op_folder <- "/Users/caitlincherryh/Documents/TestAlignmentResults/2_testFixedTrees/"
  	results_folder <- "/Users/caitlincherryh/Documents/TestAlignmentResults/2_testFixedTrees/"
  	maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/"
  	exec_folder <- "/Users/caitlincherryh/Documents/Executables/"
  	# Create a vector with all of the executable file paths
  	# To access a path: exec_paths[["name"]]
  	exec_paths <- c("3seq","iqtree","Phi","SimBac","SplitsTree.app/Contents/MacOS/JavaApplicationStub")
  	exec_paths <- paste0(exec_folder,exec_paths)
  	names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")
  	network_functions <- "code/func_split_decomposition.R"
  	run_id <- "bootstrapTest"
} else if (run_location == "nci") {
  	op_folder <- "/short/xf1/cac599/sims_output/"
  	maindir <- "/home/599/cac599/treelikeness/"
  	exec_folder <- "/home/599/cac599/treelikeness/executables/"
  	# Create a vector with all of the executable file paths
  	# To access a path: exec_paths[["name"]]
  	exec_paths <- c("3seq","iqtree","Phi","SimBac")
  	exec_paths <- paste0(exec_folder,exec_paths)
  	exec_paths <- c(exec_paths,"/home/599/cac599/splitstree4/SplitsTree")
  	names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")
  	network_functions <- "code/func_split_decomposition_nci.R"
  	run_id <- "nci2"
} else if (run_location=="soma"){
	op_folder <- "/data/caitlin/treelikeness/output/"
	results_folder <- "/data/caitlin/treelikeness/results/"
	maindir <- "/data/caitlin/treelikeness/"
	exec_paths <- c("/data/caitlin/linux_executables/3seq/3seq","/data/caitlin/linux_executables/iqtree/bin/iqtree","/data/caitlin/linux_executables/PhiPack/Phi",
	                "/data/caitlin/linux_executables/SimBac/SimBac","/data/caitlin/splitstree4/SplitsTree")
	names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")
	network_functions <- "code/func_split_decomposition.R"
	run_id <- "soma_fullSet"
}

# Set working directory
setwd(maindir)

# Source files for functions
source(paste0(maindir,network_functions))
source(paste0(maindir,"code/func_test_statistic.R"))
source(paste0(maindir,"code/func_create_alignments.R"))
source(paste0(maindir,"code/func_process_data.R"))
source(paste0(maindir,"code/func_parametric_bootstrap.R"))
tree_folder <- paste0(maindir,"trees/")

temp_time <- toc()
times <- c(times, temp_time$msg)
time_ids <- c(time_ids, (temp_time$toc - temp_time$tic)[[1]])
tic("plot1")

### Create dataframe for the final set of simulations (fixed trees)
### Each row needs to include: output_folder, n_sites, tree_age, mean_molecular_rate, sd_molecular_rate, tree1, tree2, proportion_tree2,id,rep

## For first set of plots:
# Make empty dataframe:
plot1_df <- data.frame((matrix(ncol = 8, nrow = 0)))
names(plot1_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
# Parameters that are the same for each set of trees:
output_folder <- op_folder
n_sites <- 1300
tree_age <- c(0.05, 0.1, 0.5, 1)
proportion_tree2 <- 0.5
id <- "plot1"
rep <- 1:100
tree1_vector <- c("08taxa_balanced_LHS","08taxa_balanced_LHS","08taxa_balanced_LHS","08taxa_balanced_LHS","08taxa_balanced_LHS","08taxa_balanced_LHS",
                  "08taxa_intermediate_LHS","08taxa_intermediate_LHS","08taxa_intermediate_LHS","08taxa_intermediate_LHS","08taxa_intermediate_LHS","08taxa_intermediate_LHS",
                  "08taxa_unbalanced_LHS","08taxa_unbalanced_LHS","08taxa_unbalanced_LHS","08taxa_unbalanced_LHS","08taxa_unbalanced_LHS","08taxa_unbalanced_LHS")
tree2_vector <- c("08taxa_balanced_RHS_reciprocal_close_1event","08taxa_balanced_RHS_reciprocal_divergent_1event","08taxa_balanced_RHS_reciprocal_ancient_1event",
                  "08taxa_balanced_RHS_nonreciprocal_close_1event","08taxa_balanced_RHS_nonreciprocal_divergent_1event","08taxa_balanced_RHS_nonreciprocal_ancient_1event",
                  "08taxa_intermediate_RHS_reciprocal_close_1event","08taxa_intermediate_RHS_reciprocal_divergent_1event","08taxa_intermediate_RHS_reciprocal_ancient_1event",
                  "08taxa_intermediate_RHS_nonreciprocal_close_1event","08taxa_intermediate_RHS_nonreciprocal_divergent_1event","08taxa_intermediate_RHS_nonreciprocal_ancient_1event",
                  "08taxa_unbalanced_RHS_reciprocal_close_1event","08taxa_unbalanced_RHS_reciprocal_divergent_1event","08taxa_unbalanced_RHS_reciprocal_ancient_1event",
                  "08taxa_unbalanced_RHS_nonreciprocal_close_1event","08taxa_unbalanced_RHS_nonreciprocal_divergent_1event","08taxa_unbalanced_RHS_nonreciprocal_ancient_1event")
tree_id <- 1:18
for (i in tree_id){
  tree1_temp <- tree1_vector[[i]]
  tree2_temp <- tree2_vector[[i]]
  temp_df <- expand.grid(op_folder,n_sites,tree_age,tree1_temp,tree2_temp,proportion_tree2,id,rep, stringsAsFactors = FALSE)
  names(temp_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
  plot1_df <- rbind(plot1_df,temp_df, stringsAsFactors = FALSE)
}
mclapply(1:nrow(plot1_df), phylo.fixedtrees.wrapper, plot1_df, exec_paths, tree_folder, mc.cores = 35) # mclapply for phylo with fixed trees

temp_time <- toc()
times <- c(times, temp_time$msg)
time_ids <- c(time_ids, (temp_time$toc - temp_time$tic)[[1]])
tic("plot2")

## For second set of plots:
# Make empty dataframe:
plot2_df <- data.frame((matrix(ncol = 8, nrow = 0)))
names(plot2_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
# Parameters that are the same for each set of trees:
output_folder <- op_folder
n_sites <- 1300
tree_age <- c(0.05, 0.1, 0.5, 1)
proportion_tree2 <- seq(0,0.5,0.01)
id <- "plot2"
rep <- 1:10
tree1_vector <- c("08taxa_balanced_LHS","08taxa_balanced_LHS",
                  "08taxa_intermediate_LHS","08taxa_intermediate_LHS",
                  "08taxa_unbalanced_LHS","08taxa_unbalanced_LHS")
tree2_vector <- c("08taxa_balanced_RHS_reciprocal_close_1event","08taxa_balanced_RHS_nonreciprocal_close_1event",
                  "08taxa_intermediate_RHS_reciprocal_close_1event","08taxa_intermediate_RHS_nonreciprocal_close_1event",
                  "08taxa_unbalanced_RHS_reciprocal_close_1event","08taxa_unbalanced_RHS_nonreciprocal_close_1event")
tree_id <- 1:6
for (i in tree_id){
  tree1_temp <- tree1_vector[[i]]
  tree2_temp <- tree2_vector[[i]]
  temp_df <- expand.grid(op_folder,n_sites,tree_age,tree1_temp,tree2_temp,proportion_tree2,id,rep, stringsAsFactors = FALSE)
  names(temp_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
  plot2_df <- rbind(plot2_df,temp_df, stringsAsFactors = FALSE)
}
mclapply(1:nrow(plot2_df), phylo.fixedtrees.wrapper, plot2_df, exec_paths, tree_folder, mc.cores = 35) # mclapply for phylo with fixed trees

temp_time <- toc()
times <- c(times, temp_time$msg)
time_ids <- c(time_ids, (temp_time$toc - temp_time$tic)[[1]])
tic("plot3")

## For third set of plots:
# Make empty dataframe:
plot3_df <- data.frame((matrix(ncol = 8, nrow = 0)))
names(plot3_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
# Parameters that are the same for each set of trees:
output_folder <- op_folder
n_sites <- 1300
tree_age <- c(0.05, 0.1, 0.5, 1)
proportion_tree2 <- (0.5)
id <- "plot3"
rep <- 1:100
tree1_vector <- c("32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS",
                  "32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS",
                  "32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS",
                  "32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS","32taxa_balanced_LHS")
tree2_vector <- c("32taxa_balanced_RHS_reciprocal_close_1event","32taxa_balanced_RHS_reciprocal_close_2event",
                  "32taxa_balanced_RHS_reciprocal_close_3event","32taxa_balanced_RHS_reciprocal_close_4event",
                  "32taxa_balanced_RHS_reciprocal_close_5event","32taxa_balanced_RHS_reciprocal_close_6event",
                  "32taxa_balanced_RHS_reciprocal_close_7event","32taxa_balanced_RHS_reciprocal_close_8event",
                  "32taxa_balanced_RHS_nonreciprocal_close_1event","32taxa_balanced_RHS_nonreciprocal_close_2event",
                  "32taxa_balanced_RHS_nonreciprocal_close_3event","32taxa_balanced_RHS_nonreciprocal_close_4event",
                  "32taxa_balanced_RHS_nonreciprocal_close_5event","32taxa_balanced_RHS_nonreciprocal_close_6event",
                  "32taxa_balanced_RHS_nonreciprocal_close_7event","32taxa_balanced_RHS_nonreciprocal_close_8event")
tree_id <- 1:16
for (i in tree_id){
  tree1_temp <- tree1_vector[[i]]
  tree2_temp <- tree2_vector[[i]]
  temp_df <- expand.grid(op_folder,n_sites,tree_age,tree1_temp,tree2_temp,proportion_tree2,id,rep, stringsAsFactors = FALSE)
  names(temp_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
  plot3_df <- rbind(plot3_df,temp_df, stringsAsFactors = FALSE)
}
mclapply(1:nrow(plot3_df), phylo.fixedtrees.wrapper, plot3_df, exec_paths, tree_folder, mc.cores = 35) # mclapply for phylo with fixed trees

temp_time <- toc()
times <- c(times, temp_time$msg)
time_ids <- c(time_ids, (temp_time$toc - temp_time$tic)[[1]])
tic("plot4")

## For fourth set of plots:
# Make empty dataframe:
plot4_df <- data.frame((matrix(ncol = 10, nrow = 0)))
names(plot4_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
# Parameters that are the same for each set of trees:
output_folder <- op_folder
n_sites <- 1300
tree_age <- c(0.05, 0.1, 0.5, 1)
proportion_tree2 <- seq(0,0.5,0.1)
# plot4_id <- "plot4"
plot4_id <- "TESTBOOTSTRAP"
rep <- 1:100
tree1_vector <- c("08taxa_balanced_LHS","08taxa_balanced_LHS",
                  "08taxa_intermediate_LHS","08taxa_intermediate_LHS",
                  "08taxa_unbalanced_LHS","08taxa_unbalanced_LHS")
tree2_vector <- c("08taxa_balanced_RHS_reciprocal_close_1event","08taxa_balanced_RHS_nonreciprocal_close_1event",
                  "08taxa_intermediate_RHS_reciprocal_close_1event","08taxa_intermediate_RHS_nonreciprocal_close_1event",
                  "08taxa_unbalanced_RHS_reciprocal_close_1event","08taxa_unbalanced_RHS_nonreciprocal_close_1event")
tree_id <- 1:6
for (i in tree_id){
  tree1_temp <- tree1_vector[[i]]
  tree2_temp <- tree2_vector[[i]]
  temp_df <- expand.grid(op_folder,n_sites,tree_age,tree1_temp,tree2_temp,proportion_tree2,plot4_id,rep, stringsAsFactors = FALSE)
  names(temp_df) <- c("output_folder", "n_sites", "tree_age", "tree1", "tree2", "proportion_tree2", "id", "rep")
  plot4_df <- rbind(plot4_df,temp_df, stringsAsFactors = FALSE)
}
# mclapply(1:nrow(plot4_df), phylo.fixedtrees.wrapper, plot4_df, exec_paths, tree_folder, mc.cores = 35) # mclapply for phylo with fixed trees
mclapply(1:100, phylo.fixedtrees.wrapper, plot4_df, exec_paths, tree_folder, mc.cores = 35) # mclapply for phylo with fixed trees
# Collect the folders that contain the alignments for plot4
all_folders <- list.dirs(op_folder, recursive = FALSE, full.names = TRUE) # get all the directory names in the output folder
inds <- grep(plot4_id,all_folders) # find which indexes the plot4 (bootstrap) folders are at 
plot4_folders <- all_folders[inds] # get the bootstrap folders
plot4_folders <- paste0(plot4_folders,"/") # add the slash to the end so it's a path to the directory. Bootstrap function just adds "alignment.nex" not the slash.
mclapply(plot4_folders, phylo.parametric.bootstrap, 199, exec_paths[["IQTree"]], exec_paths[["SplitsTree"]], exec_paths[["Phi"]], exec_paths[["3seq"]], mc.cores = 35) # run all the bootstraps!
mclapply(plot4_folders,phylo.collate.bootstrap, mc.cores = 35) # collate the bootstrap test statistics and calculate the p-values for the test statistics

temp_time <- toc()
times <- c(times, temp_time$msg)
time_ids <- c(time_ids, (temp_time$toc - temp_time$tic)[[1]])
tic("save parameter dataframes")

# Save the parameter dataframes
op_name <- paste0(results_folder,"plot1_input_parameters_",run_id,".csv")
write.csv(plot1_df,file=op_name)
op_name <- paste0(results_folder,"plot2_input_parameters_",run_id,".csv")
write.csv(plot2_df,file=op_name)
op_name <- paste0(results_folder,"plot3_input_parameters_",run_id,".csv")
write.csv(plot3_df,file=op_name)
op_name <- paste0(results_folder,"plot4_input_parameters_",run_id,".csv")
write.csv(plot4_df,file=op_name)

temp_time <- toc()
times <- c(times, temp_time$msg)
time_ids <- c(time_ids, (temp_time$toc - temp_time$tic)[[1]])

# save the times
time_df <- data.frame(time_ids,times)
op_name <- paste0(results_folder,"run_times_",run_id,".csv")
write.csv(time_df,file=op_name)


