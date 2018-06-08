# R code to create dataframes with a row for each simulation to run, simulate those alignments and output a csv with results of all test statistics

# Open packages
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
source(paste0(maindir,"code/split_decomposition.R"))
source(paste0(maindir,"code/test_statistic.R"))
source(paste0(maindir,"code/create_alignments.R"))

# Create a vector with all of the executable file paths
# To access a path: exec_paths[["name"]]
exec_folder <- "/Users/caitlincherryh/Documents/Executables/"
exec_paths <- c("3seq","iqtree","Phi","SimBac","SplitsTree.app/Contents/MacOS/JavaApplicationStub")
exec_paths <- paste0(exec_folder,exec_paths)
names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")

## If running 3seq remotely, need to specify the ptable before the first run so that 3seq knows where to look
#system("./3seq -f mtDNA.aln -ptable PvalueTable500 -id myFirstRun") # this command associates the Ptable with 3seq - needs a sample alignment to run correctly.

baby_phylo_df <- data.frame(matrix(nrow=0,ncol=10))
phylo_row <- c("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments2/",20,1000,0.5,1,0.1,0.1,0.5,"testBaby",1)
baby_phylo_df <- rbind(baby_phylo_df,phylo_row,phylo_row,phylo_row,stringsAsFactors=FALSE)
phylo_names <- c("output_folder","n_taxa","n_sites","birth_rate","tree_age","mean_molecular_rate","sd_molecular_rate","proportion_tree2","id","rep")
names(baby_phylo_df) <- phylo_names
baby_phylo_df[2,2]<- 180
baby_phylo_df[3,2]<- 6
baby_simbac_df <- data.frame(matrix(nrow=0,ncol=9))
simbac_row <- c("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments2/",20,1000,1000000,0.1,0,0.1,"testBaby",1)
baby_simbac_df <- rbind(baby_simbac_df,simbac_row,simbac_row,stringsAsFactors=FALSE)
simbac_names <- c("output_folder","n_taxa","n_sites","gap","internal_recombination","external_recombination","mutation_rate","id","rep")
names(baby_simbac_df) <- simbac_names

# Time to run simulations: 20 taxa, 1300 sites, all quartets = 33.883 seconds. 180 taxa, 1300 sites, all quartets = 36161.139 seconds (10 hours!). 
phylo.run1sim(baby_phylo_df[1,], exec_paths)
SimBac.run1sim(baby_simbac_df[1,],exec_paths)
#phylo.run1sim(baby_phylo_df[2,], exec_paths)
#phylo.run1sim(baby_phylo_df[3,], exec_paths)



