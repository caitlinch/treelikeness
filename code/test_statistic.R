# R program to calculate test statistics for quantifying tree likeness

# Function to call IQ-tree and return a pairwise distance matrix
iqtree.pdm <- function(iqpath,path,file){
  # Check if the tree file already exists and if it doesn't, run IQ-tree and create it
  if (file.exists(paste0(path,file,".treefile")) == FALSE){
    # Given an alignment, get a tree from IQ-tree and find the sum of the pairwise distance matrix
    system(paste0(iqpath," -s ",path,file," -nt AUTO -redo")) # call IQ-tree!
  }
  # Now the tree has been created, open it up
  tree <- read.tree(paste0(path,file,".treefile")) # read in the tree created in IQ-tree
  pdm <- cophenetic.phylo(tree) # get pairwise distance matrix for the taxa - note that this mirrors across the diagonal! (doesn't have half filled with 0's)
  pdm[upper.tri(pdm)] <- 0 # Get upper triangle and replace upper triangle coordinates (TRUE) with 0
  return(pdm)
}

# Function to open a maximum likelihood distances file from IQ tree and return a pairwise distance matrix
mldist.pdm <- function(iqpath,path,file){
  # Open the maximum likelihood distances file output when creating the tree in IQ-tree
  filename <- paste0(path,file,".mldist")
  pdm <- read.table(filename,header=FALSE,skip=1) # read in mldist file
  pdm <- pdm[2:ncol(pdm)] # remove first column (taxa names)
  pdm[upper.tri(pdm)] <- 0 # replace upper triangle coordinates (TRUE) with 0
  return(pdm)
}  

# Function to scale entries in a matrix so they sum to 1
normalise.matrix <- function(matrix){
  sum <- sum(matrix) # sum up all the pairwise distances
  matrix <- matrix/sum # scale matrix entries to add to 1
}

# Test statistic based on dividing sum of values in tree pairwise distance matrix by sum of values in alignment matrix
pdm.ratio <- function(iqpath,path,file){
  tree_pdm <- iqtree.pdm(iqpath,path,file)
  #print(tree_pdm) # to check
  tree_sum <- sum(tree_pdm) # sum up all the pairwise distances
  #print(tree_sum) # print tree sum
  
  # Calculate pairwise distances from the alignment
  alignment_pdm <- mldist.pdm(iqpath,path,file)
  #print(mldist_pdm) # to check
  alignment_sum <- sum(alignment_pdm) # sum all the pairwise distances
  #print(mldist_sum) # print tree sum
  
  # Divide sums to obtain test statistic
  ts <- tree_sum/alignment_sum
  return(ts)
}

# Test statistic based on dividing sum of values in tree pairwise distance matrix by sum of values in alignment matrix
# Adjust entries of alignment matrix and tree matrix to sum to one before dividing, divide result by number of entries in lower-triangle of matrix
normalised.pdm.ratio <- function(iqpath,path,file){
  tree_pdm <- iqtree.pdm(iqpath,path,file)
  #print(tree_pdm) # to check
  tree_pdm <- normalise.matrix(tree_pdm) # adjust matrix so all entries sum to 1
  #print(tree_pdm) # print scaled pdm matrix
  tree_sum <- sum(tree_pdm) # get new tree sum
  #print(tree_sum) # print tree sum
  
  # Calculate pairwise distances from the alignment
  # Open the maximum likelihood distances file output when creating the tree in IQ-tree
  alignment_pdm <- mldist.pdm(iqpath,path,file)
  #print(mldist_pdm) # to check
  alignment_pdm <- normalise.matrix(alignment_pdm)
  #print(mldist_pdm) # print scaled pdm matrix
  mldist_sum <- sum(mldist_pdm) # get new tree sum from normalised matrix
  #print(mldist_sum) # print tree sum
  
  # Calculate test statistic
  divisor_pdm <- (tree_pdm/mldist_pdm) # divide tree matrix by alignment matrix
  #print(divisor_pdm)
  #print(sum(divisor_pdm,na.rm = TRUE))
  ts <- sum(divisor_pdm,na.rm = TRUE) # sum all entries except the NaN ones
  fac <- factorial(length(names(tree_pdm)) - 1) # calculate the number of non-NaN entries (# of values in lower triangle)
  ts <- ts/fac # divide sum of all entries by number of entries
  return(ts)
}



# Input variables and files
alignment_path <- "/Users/caitlin/Repositories/treelikeness/raw_data" # folder where alignment is located
alignment_paths <- list.dirs(alignment_path)
# alignment_paths <- paste0(alignment_paths[2:length(alignment_paths)],"/") # to run all alignments in directory
alignment_paths <- paste0(alignment_paths[2:11],"/") # to run 10 alignments
alignment_file <- "alignment.nex" # name of alignment 
iqtree_path <- "/Applications/iqtree/bin/iqtree" # location of IQ-tree program

# Open packages
library(ape)

# Run test statistic(s) 
names <- c()
all_ts <- c()
all_ts_norm <- c()
for (alignment in alignment_paths){
  print(alignment)
  names <- c(names,alignment)
  ts1 <- pdm.ratio(iqtree_path,alignment,alignment_file)
  print(ts1)
  ts2 <- normalised.pdm.ratio(iqtree_path,alignment,alignment_file)
  print(ts2)
  all_ts <- c(all_ts,ts1)
  all_ts_norm <- c(all_ts_norm,ts2)
}
  

