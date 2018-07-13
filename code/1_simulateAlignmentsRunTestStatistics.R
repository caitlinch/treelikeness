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

run_location <- "nci"
run_location <- "mac"

if (run_location == "mac"){
  op_folder <- "/Users/caitlincherryh/Documents/TestAlignmentResults/1_mainrun/"
  maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/"
  exec_folder <- "/Users/caitlincherryh/Documents/Executables/"
  network_functions <- "code/func_split_decomposition.R"
} else if (run_location == "nci") {
  op_folder <- "/short/xf1/cac599/sims_output/"
  maindir <- "/home/599/cac599/treelikeness/"
  exec_folder <- "/home/599/cac599/treelikeness/executables/"
  network_functions <- "code/func_split_decomposition_nci.R"
}

# Set working directory
setwd(maindir)

# Source files for functions
source(paste0(maindir,network_functions))
source(paste0(maindir,"code/func_test_statistic.R"))
source(paste0(maindir,"code/func_create_alignments.R"))
source(paste0(maindir,"code/func_process_data.R"))

# Create a vector with all of the executable file paths
# To access a path: exec_paths[["name"]]
exec_paths <- c("3seq","iqtree","Phi","SimBac","SplitsTree.app/Contents/MacOS/JavaApplicationStub")
exec_paths <- paste0(exec_folder,exec_paths)
names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")

## If running 3seq remotely, need to specify the ptable before the first run so that 3seq knows where to look
system(exec_paths[["3seq"]],"-f mtDNA.aln -ptable PvalueTable500 -id myFirstRun") # this command associates the Ptable with 3seq - needs a sample alignment to run correctly.

# Columns required for the rows of parameters to run the simulations
phylo_df_names <- c("output_folder","n_taxa","n_sites","birth_rate","tree_age","mean_molecular_rate","sd_molecular_rate","proportion_tree2","id","rep")
simbac_df_names <- c("output_folder","n_taxa","n_sites","gap","internal_recombination","external_recombination","mutation_rate","id","rep")

# lapply(1:nrow(external_df),SimBac.rowWrapper,dataframe = external_df, program_paths = exec_paths) # lapply for SimBac
# lapply(1:nrow(phylo_df),phylo.rowWrapper,dataframe = phylo_df, program_paths = exec_paths) # lapply for phylo

# For external recombination
# Create dataframe
output_folder <- c(op_folder)
n_taxa <- c(5, 10, 20, 40, 80, 160)
n_sites <- c(1300)
gap <- c(1000000)
internal_recombination <- c(0)
external_recombination <- c(0.000000001, 0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1)
mutation_rate <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
id <- c("external")
rep <- c(1:10)
external_df <- expand.grid(output_folder,n_taxa,n_sites,gap,internal_recombination,external_recombination,mutation_rate,id,rep, stringsAsFactors = FALSE)
names(external_df) <- c("output_folder","n_taxa","n_sites","gap","internal_recombination","external_recombination","mutation_rate","id","rep")
# run simulations
#lapply(1:nrow(internal_df),SimBac.rowWrapper,dataframe = external_df, program_paths = exec_paths)
lapply(1:2,SimBac.rowWrapper,dataframe = external_df, program_paths = exec_paths)

# For internal recombination
# Create dataframe
output_folder <- c(op_folder)
n_taxa <- c(5, 10, 20, 40, 80, 160)
n_sites <- c(1300)
gap <- c(1000000)
internal_recombination <- c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)
external_recombination <- c(0)
mutation_rate <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
id <- c("internal")
rep <- c(1:10)
internal_df <- expand.grid(output_folder,n_taxa,n_sites,gap,internal_recombination,external_recombination,mutation_rate,id,rep, stringsAsFactors = FALSE)
names(internal_df) <- c("output_folder","n_taxa","n_sites","gap","internal_recombination","external_recombination","mutation_rate","id","rep")
# run simulations
#lapply(1:nrow(internal_df),SimBac.rowWrapper,dataframe = internal_df, program_paths = exec_paths)
  
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
#lapply(1:nrow(phylo_df),phylo.rowWrapper,dataframe = phylo_df, program_paths = exec_paths)



