# R functions to calculate 3 test statistics for quantifying tree likeness

# Function to call IQ-tree and run it
call.IQTREE <- function(iqtree_path,alignment_path){
  # Check if the tree file already exists and if it doesn't, run IQ-tree and create it
  if (file.exists(paste0(alignment_path,".treefile")) == FALSE){
    # Given an alignment, get a tree from IQ-tree and find the sum of the pairwise distance matrix
    system(paste0(iqtree_path," -s ",alignment_path," -nt AUTO -redo")) # call IQ-tree!
  }
}

# Function get a tree from IQ-tree and return a pairwise distance matrix
iqtree.pdm <- function(iqpath,path){
  # Check if the tree file already exists and if it doesn't, run IQ-tree and create it
  call.IQTREE(iqpath,path)
  # Now the tree has been created, open it up
  tree <- read.tree(paste0(path,".treefile")) # read in the tree created in IQ-tree
  pdm <- cophenetic.phylo(tree) # get pairwise distance matrix for the taxa - note that this mirrors across the diagonal! (doesn't have half filled with 0's)
  pdm[upper.tri(pdm)] <- 0 # Get upper triangle and replace upper triangle coordinates (TRUE) with 0
  return(pdm)
}

# Function to open a maximum likelihood distances file from IQ tree and return a pairwise distance matrix
mldist.pdm <- function(iqpath,path){
  # Open the maximum likelihood distances file output when creating the tree in IQ-tree
  filename <- paste0(path,".mldist")
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

# Test statistic 1: based on dividing sum of values in tree pairwise distance matrix by sum of values in alignment matrix
pdm.ratio <- function(iqpath,path){
  tree_pdm <- iqtree.pdm(iqpath,path)
  #print(tree_pdm) # to check
  tree_sum <- sum(tree_pdm) # sum up all the pairwise distances
  #print(tree_sum) # print tree sum
  
  # Calculate pairwise distances from the alignment
  alignment_pdm <- mldist.pdm(iqpath,path)
  #print(mldist_pdm) # to check
  alignment_sum <- sum(alignment_pdm) # sum all the pairwise distances
  #print(mldist_sum) # print tree sum
  
  # Divide sums to obtain test statistic
  ts <- tree_sum/alignment_sum
  return(ts)
}

# Test statistic 2: sum of the absolute differences of the normalised distance matrices
# Adjust entries of alignment matrix and tree matrix to sum to one before dividing, divide result by number of entries in lower-triangle of matrix
normalised.pdm.diff.sum <- function(iqpath,path){
  # Calculate the pairwise differences matrix for the tree and normalise it
  tree_pdm <- iqtree.pdm(iqpath,path)
  tree_pdm <- normalise.matrix(tree_pdm) # adjust matrix so all entries sum to 1
  
  # Calculate pairwise distances from the alignment and normalise it
  # Open the maximum likelihood distances file output when creating the tree in IQ-tree
  alignment_pdm <- mldist.pdm(iqpath,path)
  alignment_pdm <- normalise.matrix(alignment_pdm)
  
  # Find the absolute difference between the tree pdm and the alignment pdm
  # gives the relative proportion that each element takes up
  diff_pdm <- abs(tree_pdm - alignment_pdm)
  # sum the difference matrix to get the test statistic
  ts <- sum(diff_pdm)
  # return the test statistic
  return(ts)
}

split.decomposition.statistic <- function(taxa_names,distance_matrix,phylogenetic_tree){
  # Open/get matrix
  # Get split decomposition
  # Get tree
  # Which splits in the decomposition are in the tree (e.g. which are monophyletic)?
  # Sum split weights of splits in tree
  # Divide by sum of all split weights
  # Return test statistic
}

