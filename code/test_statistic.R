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
mldist.pdm <- function(path){
  # Open the maximum likelihood distances file output from creating the tree in IQ-tree
  filename <- paste0(path,".mldist")
  pdm <- read.table(filename,header=FALSE,skip=1) # read in mldist file
  pdm <- pdm[2:ncol(pdm)] # remove first column (taxa names)
  pdm[upper.tri(pdm)] <- 0 # replace upper triangle coordinates (TRUE) with 0
  return(pdm)
}  

# Function to open a maximum likelihood distances file from IQ tree and return the taxa names
mldist.taxa <- function(path){
  # Open the maximum likelihood distances file output from creating the tree in IQ-tree
  filename <- paste0(path,".mldist")
  pdm <- read.table(filename,header=FALSE,skip=1) # read in mldist file
  taxa_names <- as.character(pdm[,1]) # extract first column (list of taxa names)
  return(taxa_names)
}

# Function to scale entries in a matrix so they sum to 1
normalise.matrix <- function(matrix){
  sum <- sum(matrix) # sum up all the pairwise distances
  matrix <- matrix/sum # scale matrix entries to add to 1
}

# Function to open an IQ-TREE generated tree given the name of the alignment
open.tree <- function(path){
  tree_file <- paste0(path,".treefile")
  tree <- read.tree(tree_file)
  return(tree)
}

# Test statistic 1: based on dividing sum of values in tree pairwise distance matrix by sum of values in alignment matrix
pdm.ratio <- function(iqpath,path){
  tree_pdm <- iqtree.pdm(iqpath,path)
  tree_sum <- sum(tree_pdm) # sum up all the pairwise distances
  
  # Calculate pairwise distances from the alignment
  alignment_pdm <- mldist.pdm(path)
  alignment_sum <- sum(alignment_pdm) # sum all the pairwise distances
  
  # Divide sums to obtain test statistic
  ts <- tree_sum/alignment_sum
  return(ts)
}

# Test statistic 2a: sum of the absolute differences of the normalised distance matrices
# Adjust entries of alignment matrix and tree matrix to sum to one before dividing, divide result by number of entries in lower-triangle of matrix
normalised.pdm.diff.sum <- function(iqpath,path){
  # Calculate the pairwise differences matrix for the tree and normalise it
  tree_pdm <- iqtree.pdm(iqpath,path)
  tree_pdm <- normalise.matrix(tree_pdm) # adjust matrix so all entries sum to 1
  
  # Calculate pairwise distances from the alignment and normalise it
  # Open the maximum likelihood distances file output when creating the tree in IQ-tree
  alignment_pdm <- mldist.pdm(path)
  alignment_pdm <- normalise.matrix(alignment_pdm)
  
  # Find the absolute difference between the tree pdm and the alignment pdm
  # gives the relative proportion that each element takes up
  diff_pdm <- abs(tree_pdm - alignment_pdm)
  # sum the difference matrix to get the test statistic
  ts <- sum(diff_pdm)
  # return the test statistic
  return(ts)
}

# Test statistic 2b: mean of the absolute differences of the normalised distance matrices
# Adjust entries of alignment matrix and tree matrix to sum to one before dividing, divide result by number of entries in lower-triangle of matrix
# More comparable then test statistic 2a as datasets with different numbers of taxa and therefore different matrix sizes will be normalised
#     to the number of taxa (because the mean is just the sum divided by the number of values).
normalised.pdm.diff.mean <- function(iqpath,path){
  # Calculate the pairwise differences matrix for the tree and normalise it
  tree_pdm <- iqtree.pdm(iqpath,path)
  tree_pdm <- normalise.matrix(tree_pdm) # adjust matrix so all entries sum to 1
  
  # Calculate pairwise distances from the alignment and normalise it
  # Open the maximum likelihood distances file output when creating the tree in IQ-tree
  alignment_pdm <- mldist.pdm(path)
  alignment_pdm <- normalise.matrix(alignment_pdm)
  
  # Find the absolute difference between the tree pdm and the alignment pdm
  # gives the relative proportion that each element takes up
  diff_pdm <- abs(tree_pdm - alignment_pdm)
  # sum the difference matrix to get the test statistic
  ts <- mean(diff_pdm)
  # return the test statistic
  return(ts)
}

# Test statistic 3: proportion of all split weights present in the tree
# Find which splits are in the tree and sum those split weights, divide by sum of all split weights
split.decomposition.manual.statistic <- function(iq_path,path){
  ## Run IQ-tree if it hasn't already been run
  call.IQTREE(iqpath,path) # path refers to an alignments
  
  ## Calculate the split decomposition
  # Open pairwise distance matrix from IQ-TREE (use mldist matrix because it uses a model to compensate for saturation)
  pdm <- mldist.pdm(path)
  taxa <- mldist.taxa(path)
  # Get split decomposition
  splits <- split.decomposition(taxa, pdm, threshold = 1)
 
  ## Open the tree estimated by IQ-TREE
  tree <- open.tree(path)
  
  ## Identify which splits in the split decomposition are in the tree (e.g. which are monophyletic)
  # Initialise sum of weights-of-splits-in-tree
  intree_sum <- 0
  # Initialise sum of weights of all splits
  all_sum <- 0
  # Iterate through each split
  splits_list <- splits[[1]]
  for (s in splits_list){
    # Extract the subsets of taxa for the split
    ss1 <- s$partition$a
    ss2 <- s$partition$b
    # Test whether the split is in the tree
    ss1_mono <- is.monophyletic(tree,ss1)
    ss2_mono <- is.monophyletic(tree,ss2)
    # If both are true then add the split weight to the sum of split weights from splits in the tree
    # (if one is true, the other must be true - reduces to (A,B) where A and B are the two monophyletic clades)
    # Checks both just to be sure 
    if (ss1_mono == TRUE & ss2_mono == TRUE){
      intree_sum <- intree_sum + s$isolation_index
    }
    # Add the split weight to the all_sum regardless of whether the split is in the tree
    # to have running total of the sum of split weights
    all_sum <- all_sum + s$isolation_index
  }
  
  ## Calculate test statistic
  # Divid sum of weights-of-splits-in-tree by sum of weights-of-all-splits
  ts <- intree_sum/all_sum
  # Return test statistic
  return(ts)
}

# Test statistic 3: proportion of all split weights present in the tree
# Find which splits are in the tree and sum those split weights, divide by sum of all split weights#
# network_algorithm - either "split decomposition" or "neighbournet" - defines which transformation will be applied
#       to the alignment to turn it into a network in SplitsTree
SplitsTree.decomposition.statistic <- function(iqpath, splitstree_path, path,network_algorithm){
  # Run IQ-tree if it hasn't already been run
  call.IQTREE(iqpath,path) # path = path to alignment
  # Calculate the split decomposition
  call.SplitsTree(splitstree_path,path,network_algorithm)
  # Retrieve the file name for the splits output file
  splits.filepath <- splits.filename(path)
  # Extract the splits 
  splits <- read.nexus.splits(splits.filepath) # Open the splits from SplitsTree
  ## Open the tree estimated by IQ-TREE
  tree <- open.tree(path)
  
  # Test each split to see whether it's in the tree
  # Make sure to feed in each line using [x] into the apply function, not [[x]] or the attributes (taxa names, weights) won't be passed to the test.monophyly function
  tree_ii_sum <- 0 # create a vector to store split weights
  for (i in 1:length(splits)){
    # Iterate through each of the rows in the splits dataframe
    test_ii <- test.monophyly(splits[i],tree) # test for monophyly
    tree_ii_sum <- tree_ii_sum + test_ii # add weight to running sum (or 0 if split not in tree)
  }
  # Get the sum of all split weights
  all_ii_sum <- sum.all.ii(splits)
  # Divide the sum of split weights in the tree by the sum of all split weights
  ts <- tree_ii_sum / all_ii_sum
  # Return the test statistic result
  return(ts)
}

# Function to test monophyly of a split, and if it's monophyletic, return the isolation index.
test.monophyly <- function(split, tree){
  # Extract the bipartition subsets from the tree
  ss1 <- split[[1]] # get the indices of all taxa in the split
  taxa <- attr(split,"labels") # get the names of all taxa in the tree
  ss1_taxa <- taxa[ss1] # use the split indices to get the taxa names from the split
  ss2_taxa <- setdiff(taxa,ss1_taxa) # take the elements from taxa not in subset 1 to get the elements in ss2
  # Test whether the split is in the tree by seeing whether ss1 and ss2 are monophyletic
  ss1_mono <- is.monophyletic(tree,ss1_taxa)
  ss2_mono <- is.monophyletic(tree,ss2_taxa)
  # If both are true then add the split weight to the sum of split weights from splits in the tree
  # (if one is true, the other must be true - reduces to (A,B) where A and B are the two monophyletic clades)
  # Checks both just to be sure 
  # If they're not both tre
  if (ss1_mono == TRUE & ss2_mono == TRUE){
    ii <- attr(split, "weights") # set the isolation index to the split weight
  } else {
    ii <- 0 # set the isolation index to 0 as this split is not in the tree
  }
  return(ii)
}

# Function to sum all weights from a splits nexus file
sum.all.ii <- function(splits){
  total_ii <- sum(attr(splits,"weights"))
  return(total_ii)
}
