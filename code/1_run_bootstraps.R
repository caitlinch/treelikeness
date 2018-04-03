# Code to run the parametric bootstraps for the relevant datasets

# Open packages
library(ape)
library(TreeSim)
library(phytools)
library(phangorn)
library(readtext) # for reading iq-tree output

# Set working directory
maindir <- "/Users/caitlin/Repositories/treelikeness/"
setwd(maindir)

# Source files for functions
source(paste0(maindir,"code/split_decomposition.R"))
source(paste0(maindir,"code/parametric_bootstrap.R"))
source(paste0(maindir,"code/test_statistic.R"))


## Test code for split decomposition functions
# # Enter distance matrix and taxa labels for practicing
# d <- t(matrix(c(0,0,0,0,0,0,0,4.0,0,0,0,0,0,0,5.0,1.0,0,0,0,0,0,7.0,3.0,2.0,0,0,0,0,13.0,9.0,8.0,6.0,0,0,
#                 0,8.0,12.0,13.0,11.0,5.0,0,0,6.0,10.0,11.0,13.0,7.0,2.0,0),
#               nrow = 7, ncol = 7)) # transpose as R filles by columns first not by rows first
# taxa <- c("A","B","C","D","E","F","G")
# rownames(d) <- taxa # label rownames with taxa
# colnames(d) <- taxa # label colnames with taxa
# a = split_decomposition(taxa,d)



## Test code to call test statistics and run them
# # Input variables and files
alignment_path <- "/Users/caitlin/Repositories/treelikeness/raw_data" # folder where alignment is located
alignment_paths <- list.dirs(alignment_path)
alignment_paths <- paste0(alignment_paths[2:length(alignment_paths)],"/") # to run all alignments in directory
alignment_file <- "alignment.nex" # name of alignment 
iqtree_path <- "/Applications/iqtree/bin/iqtree" # location of IQ-tree program

for (alignment in alignment_paths){
   print(alignment)
   system(paste0(iqtree_path," -s ",alignment,alignment_file," -nt AUTO -redo"))
}

### User Input Parameters
nbootstrap <- 3 # the number of parametric bootstraps you'd like to run
filepath <- "" # The folder containing all of the datasets of interest


### To run the parametric bootstrap 99/999 times
# generate a list of ids for each simulation
ids <- sprintf("%04d",1:nbootstrap)
for (id in ids){
  # call the parametric bootstrap function on the folder of interest
}









