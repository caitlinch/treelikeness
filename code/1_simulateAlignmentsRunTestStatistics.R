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
# library(tictoc) # library for measuring timings! tic("label") to start, toc() to stop

# Remember to have downloaded and tested the following programs: SplitsTree4, SimBac, IQ-Tree, PhiPack and 3Seq. 
# 3Seq must have been associated with a P-value table for it to run properly
# If it's not working on soma, you (Caitlin) might have forgotten about the LD_LIBRARY_PATH again: run these two lines of code
# LD_LIBRARY_PATH=$LD_LIBRARY_PATH:usr/local/lib:usr/lib/x86_64-linux-gnu
# export LD_LIBRARY_PATH

# run_location <- "nci"
# run_location <- "mac"
run_location <- "soma"
run_id <- "soma1"

if (run_location == "mac"){
  	op_folder <- "/Users/caitlincherryh/Documents/TestAlignmentResults/2_testFixedTrees/"
  	maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/"
  	exec_folder <- "/Users/caitlincherryh/Documents/Executables/"
  	# Create a vector with all of the executable file paths
  	# To access a path: exec_paths[["name"]]
  	exec_paths <- c("3seq","iqtree","Phi","SimBac","SplitsTree.app/Contents/MacOS/JavaApplicationStub")
  	exec_paths <- paste0(exec_folder,exec_paths)
  	names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")
  	network_functions <- "code/func_split_decomposition.R"
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
} else if (run_location=="soma"){
	op_folder <- "/data/caitlin/treelikeness/output/"
	results_folder <- "/data/caitlin/treelikeness/results/"
	maindir <- "/data/caitlin/treelikeness/"
	exec_paths <- c("/data/caitlin/linux_executables/3seq/3seq","/data/caitlin/linux_executables/iqtree/bin/iqtree","/data/caitlin/linux_executables/PhiPack/Phi",
	                "/data/caitlin/linux_executables/SimBac/SimBac","/data/caitlin/splitstree4/SplitsTree")
	names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")
	network_functions <- "code/func_split_decomposition.R"
}

# Set working directory
setwd(maindir)

# Source files for functions
source(paste0(maindir,network_functions))
source(paste0(maindir,"code/func_test_statistic.R"))
source(paste0(maindir,"code/func_create_alignments.R"))
source(paste0(maindir,"code/func_process_data.R"))
tree_folder <- paste0(maindir,"trees/")

# Columns required for the rows of parameters to run the simulations
phylo_df_names <- c("output_folder","n_taxa","n_sites","birth_rate","tree_age","mean_molecular_rate","sd_molecular_rate","proportion_tree2","id","rep")
simbac_df_names <- c("output_folder","n_taxa","n_sites","gap","internal_recombination","external_recombination","mutation_rate","id","rep")

# lapply(1:nrow(external_df),SimBac.rowWrapper,dataframe = external_df, program_paths = exec_paths) # lapply for SimBac
# lapply(1:nrow(phylo_df),phylo.rowWrapper,dataframe = phylo_df, program_paths = exec_paths) # lapply for phylo

# For external recombination
# Create dataframe
output_folder <- c(op_folder)
n_taxa <- c(5, 10, 20, 40, 80, 160) # removed 5
n_sites <- c(1300)
gap <- c(1000000)
internal_recombination <- c(0)
external_recombination <- c(0.000000001, 0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1) #removed 0.000000001, 0.00000001, 0.0000001, 0.000001, 0.00001
mutation_rate <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05) # removed 0.0001, 0.0005
id <- c("external")
rep <- c(1:10)
external_df <- expand.grid(output_folder,n_taxa,n_sites,gap,internal_recombination,external_recombination,mutation_rate,id,rep, stringsAsFactors = FALSE)
names(external_df) <- c("output_folder","n_taxa","n_sites","gap","internal_recombination","external_recombination","mutation_rate","id","rep")
# run simulations
mclapply(1:nrow(external_df),SimBac.rowWrapper,dataframe = external_df, program_paths = exec_paths, mc.cores = 35)
#mclapply(1:2,SimBac.rowWrapper,dataframe = external_df, program_paths = exec_paths, mc.cores = 35)

# For internal recombination
# Create dataframe
output_folder <- c(op_folder)
n_taxa <- c(5, 10, 20, 40, 80, 160) #removed 5
n_sites <- c(1300)
gap <- c(1000000)
internal_recombination <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
external_recombination <- c(0)
mutation_rate <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05) # removed 0.0001, 0.0005
id <- c("internal")
rep <- c(1:10)
internal_df <- expand.grid(output_folder,n_taxa,n_sites,gap,internal_recombination,external_recombination,mutation_rate,id,rep, stringsAsFactors = FALSE)
names(internal_df) <- c("output_folder","n_taxa","n_sites","gap","internal_recombination","external_recombination","mutation_rate","id","rep")
# run simulations
mclapply(1:nrow(internal_df),SimBac.rowWrapper,dataframe = internal_df, program_paths = exec_paths, mc.cores = 35)
  
# Create phylogenetic sims dataframe
output_folder <- c(op_folder)
n_taxa <- c(5, 10, 20, 40, 80, 160)
n_sites <- c(1300)
birth_rate <- c(0.5)
tree_age <- c(1)
mean_molecular_rate <- c(1, 0.5, 0.1, 0.05)
sd_molecular_rate <- c(0.1)
proportion_tree2 <- c(seq(0,0.5,0.01))
id <- c("2trees") # can't use phylo as an ID as the word is included in title of all sims made using the phylogenetic approach
rep <- c(1:10)
phylo_df <- expand.grid(output_folder,n_taxa,n_sites,birth_rate,tree_age,mean_molecular_rate,sd_molecular_rate,proportion_tree2,id,rep, stringsAsFactors = FALSE)
names(phylo_df) <- c("output_folder","n_taxa","n_sites","birth_rate","tree_age","mean_molecular_rate","sd_molecular_rate","proportion_tree2","id","rep")
# Run the simulations
mclapply(1:nrow(phylo_df),phylo.rowWrapper,dataframe = phylo_df, program_paths = exec_paths, mc.cores = 35)

# Save the parameter dataframes
op_name <- paste0(results_folder,"external_input_parameters_",run_id,".csv")
write.csv(external_df,file=op_name)
op_name <- paste0(results_folder,"internal_input_parameters_",run_id,".csv")
write.csv(internal_df,file=op_name)
op_name <- paste0(results_folder,"2trees_input_parameters_",run_id,".csv")
write.csv(phylo_df,file=op_name)


# test for fixed trees
row <- c(op_folder,1300,1,0.1,0.1,"16taxa_intermediate_LHS","16taxa_intermediate_RHS_nonreciprocal_divergent_2event",0.5,"test",1)
names(row) <- c("output_folder","n_sites","tree_age","mean_molecular_rate","sd_molecular_rate","tree1","tree2","proportion_tree2","id","rep")
