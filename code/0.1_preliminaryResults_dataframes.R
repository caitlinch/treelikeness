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

# Set working directory
maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # for work computer
setwd(maindir)

# Source files for functions
source(paste0(maindir,"code/func_split_decomposition.R"))
source(paste0(maindir,"code/func_test_statistic.R"))
source(paste0(maindir,"code/func_create_alignments.R"))
source(paste0(maindir,"code/func_process_data.R"))

# Create a vector with all of the executable file paths
# To access a path: exec_paths[["name"]]
exec_folder <- "/Users/caitlincherryh/Documents/Executables/"
exec_paths <- c("3seq","iqtree","Phi","SimBac","SplitsTree.app/Contents/MacOS/JavaApplicationStub")
exec_paths <- paste0(exec_folder,exec_paths)
names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")

## If running 3seq remotely, need to specify the ptable before the first run so that 3seq knows where to look
#system("./3seq -f mtDNA.aln -ptable PvalueTable500 -id myFirstRun") # this command associates the Ptable with 3seq - needs a sample alignment to run correctly.

# Columns required for the rows of parameters to run the simulations
phylo_df_names <- c("output_folder","n_taxa","n_sites","birth_rate","tree_age","mean_molecular_rate","sd_molecular_rate","proportion_tree2","id","rep")
simbac_df_names <- c("output_folder","n_taxa","n_sites","gap","internal_recombination","external_recombination","mutation_rate","id","rep")

# lapply(1:nrow(external_df),SimBac.rowWrapper,dataframe = external_df, program_paths = exec_paths) # lapply for SimBac
# lapply(1:nrow(phylo_df),phylo.rowWrapper,dataframe = phylo_df, program_paths = exec_paths) # lapply for phylo

# Build some test dataframes
output_folder <- rep("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments3/",20)
n_taxa <- rep(20,20)
n_sites <- rep(1300,20)
gap <- rep(1000000,20)
internal_recombination <- rep(c(0.00,0.00,0.00,0.00),5)
external_recombination <- rep(c(0.00,0.01,0.05,0.10),5)
mutation_rate <- rep(0.01,20)
id <- rep(c("external"),20)
rep <- c(rep(1,4), rep(2,4),rep(3,4),rep(4,4),rep(5,4))
external_df <- data.frame(output_folder,n_taxa,n_sites,gap,internal_recombination,external_recombination,mutation_rate,id,rep, stringsAsFactors = FALSE)
lapply(1:nrow(external_df),SimBac.rowWrapper,dataframe = external_df, program_paths = exec_paths)

output_folder <- rep("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments3/",20)
n_taxa <- rep(20,20)
n_sites <- rep(1300,20)
gap <- rep(1000000,20)
internal_recombination <- rep(c(0.00, 0.01, 0.10,0.20),5)
external_recombination <- rep(c(0.00,0.00,0.00,0.0),5)
mutation_rate <- rep(0.01,20)
id <- rep(c("internal"),20)
rep <- c(rep(1,4), rep(2,4),rep(3,4),rep(4,4),rep(5,4))
internal_df <- data.frame(output_folder,n_taxa,n_sites,gap,internal_recombination,external_recombination,mutation_rate,id,rep, stringsAsFactors = FALSE)
lapply(1:nrow(internal_df),SimBac.rowWrapper,dataframe = internal_df, program_paths = exec_paths)

output_folder <- rep("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments3/",20)
n_taxa <- rep(20,20)
n_sites <- rep(1300,20)
gap <- rep(1000000,20)
internal_recombination <- rep(c(0.00, 0.01, 0.10,0.20),5)
external_recombination <- rep(c(0.00,0.00,0.00,0.0),5)
mutation_rate <- rep(0.01,20)
id <- rep(c("mutation"),20)
rep <- c(rep(1,4), rep(2,4),rep(3,4),rep(4,4),rep(5,4))
mutation_df <- data.frame(output_folder,n_taxa,n_sites,gap,internal_recombination,external_recombination,mutation_rate,id,rep, stringsAsFactors = FALSE)
lapply(1:nrow(mutation_df),SimBac.rowWrapper,dataframe = mutation_df, program_paths = exec_paths)

output_folder <- rep("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments3/",255)
n_taxa <- rep(20,255)
n_sites <- rep(1300,255)
birth_rate <- rep(0.5,255)
tree_age <- rep(1,255)
mean_molecular_rate <- rep(0.1,255)
sd_molecular_rate <- rep(0.1,255)
proportion_tree2 <- rep(seq(0,0.5,0.01),5)
id <- rep(c("2trees"),255) # can't use phylo as an ID as the word is included in title of all sims made using the phylogenetic approach
rep <- c(rep(1,51), rep(2,51),rep(3,51),rep(4,51),rep(5,51))
phylo_df <- data.frame(output_folder,n_taxa,n_sites,birth_rate,tree_age,mean_molecular_rate,sd_molecular_rate,proportion_tree2,id,rep, stringsAsFactors = FALSE)
lapply(1:nrow(phylo_df),phylo.rowWrapper,dataframe = phylo_df, program_paths = exec_paths)

output_folder <- rep("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments3/",204)
n_taxa <- rep(20,204)
n_sites <- rep(1300,204)
birth_rate <- rep(0.5,204)
tree_age <- rep(1,204)
mean_molecular_rate <- rep(c(0.05,0.1,0.5,1),51)
sd_molecular_rate <- rep(0.1,204)
proportion_tree2 <- rep(seq(0,0.5,0.01),4)
id <- rep(c("molRate"),204)
rep <- c(rep(1,51), rep(2,51),rep(3,51),rep(4,51))
depth_df <- data.frame(output_folder,n_taxa,n_sites,birth_rate,tree_age,mean_molecular_rate,sd_molecular_rate,proportion_tree2,id,rep, stringsAsFactors = FALSE)
lapply(1:nrow(depth_df),phylo.rowWrapper,dataframe = depth_df, program_paths = exec_paths)

output_folder <- rep("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments3/",1)
n_taxa <- rep(6,1)
n_sites <- rep(1300,1)
birth_rate <- rep(0.5,1)
tree_age <- rep(1,1)
mean_molecular_rate <- rep(0.1,1)
sd_molecular_rate <- rep(0.1,1)
proportion_tree2 <- rep(0.25,1)
id <- rep(c("CheckTestStats"),1)
rep <- rep(1,1)
check_df <- data.frame(output_folder,n_taxa,n_sites,birth_rate,tree_age,mean_molecular_rate,sd_molecular_rate,proportion_tree2,id,rep, stringsAsFactors = FALSE)
lapply(1:nrow(check_df),phylo.rowWrapper,dataframe = check_df, program_paths = exec_paths)