# R functions to calculate 3 test statistics for quantifying tree likeness

library(TreeSim)
library(phytools)
library(seqinr)
library(ape)
library(phangorn)

# Function to call IQ-tree and run it
call.IQTREE <- function(iqtree_path,alignment_path){
  # Check if the tree file already exists and if it doesn't, run IQ-tree and create it
  if (file.exists(paste0(alignment_path,".treefile")) == FALSE){
    # Given an alignment, get a tree from IQ-tree and find the sum of the pairwise distance matrix
    system(paste0(iqtree_path," -s ",alignment_path," -nt 1 -redo -safe")) # call IQ-tree!
  }
}

call.IQTREE.quartet <- function(iqtree_path,alignment_path,nsequences){
  # For this alignment, check if the IQ-Tree log file, the treefile OR the likelihood map for this alignment exist
  # If any of those files don't exist, run IQ-Tree
  # This way if IQ-Tree failed to run properly, or if it ran before without doing the likelihood mapping, it will rerun here
  if (file.exists(paste0(alignment_path,".iqtree")) == FALSE){
    # Given an alignment, get a tree from IQ-tree and find the sum of the pairwise distance matrix
    # Specify -lmap with 25 times the number of sequences, so that each sequence is covered ~100 times in the quartet sampling
    nquartet <- 25*as.numeric(nsequences)
    system(paste0(iqtree_path," -s ",alignment_path," -nt 1 -lmap ",nquartet," -redo -safe")) # call IQ-tree!
  } else if (file.exists(paste0(alignment_path,".lmap.eps")) == FALSE){
    # Given an alignment, get a tree from IQ-tree and find the sum of the pairwise distance matrix
    # Specify -lmap with 25 times the number of sequences, so that each sequence is covered ~100 times in the quartet sampling
    nquartet <- 25*as.numeric(nsequences)
    system(paste0(iqtree_path," -s ",alignment_path," -nt 1 -lmap ",nquartet," -redo -safe")) # call IQ-tree!
  } else if (file.exists(paste0(alignment_path,".treefile")) == FALSE){
    # Given an alignment, get a tree from IQ-tree and find the sum of the pairwise distance matrix
    # Specify -lmap with 25 times the number of sequences, so that each sequence is covered ~100 times in the quartet sampling
    nquartet <- 25*as.numeric(nsequences)
    system(paste0(iqtree_path," -s ",alignment_path," -nt 1 -lmap ",nquartet," -redo -safe")) # call IQ-tree!
  } 
}

call.IQTREE.quartet.bootstrap <- function(iqtree_path,alignment_path,nsequences){
  # Check if the tree file already exists and if it doesn't, run IQ-tree and create it
  if (file.exists(paste0(alignment_path,".treefile")) == FALSE){
    # Given an alignment, get a tree from IQ-tree and find the sum of the pairwise distance matrix
    # Specify -lmap with 25 times the number of sequences, so that each sequence is covered ~100 times in the quartet sampling
    nquartet <- 25*as.numeric(nsequences)
    system(paste0(iqtree_path," -s ",alignment_path," -nt 1 -lmap ",nquartet," -redo -m JC -safe")) # call IQ-tree!
  }
}

# Function get a tree from IQ-tree and return a pairwise distance matrix
iqtree.pdm <- function(iqpath,path){
  # Check if the tree file already exists and if it doesn't, run IQ-tree and create it
  call.IQTREE(iqpath,path)
  # Now the tree has been created, open it up
  tree <- read.tree(paste0(path,".treefile")) # read in the tree created in IQ-tree
  pdm <- cophenetic.phylo(tree) # get pairwise distance matrix for the taxa - note that this mirrors across the diagonal! (doesn't have half filled with 0's)
  # After using cophenetic.phylo on the tree from IQ-Tree, need to rearrange the matrix
  # So that when doing the subtractions, the pairwise distances will be coming from the same location (eg AB and AB, not BF and EH)
  # matrix is not in a helpful order (orders based on taxa position in tree, but I want numerical)
  # to sort cophenetic.phylo (matrix from the tree)
  first_taxa_name <- rownames(pdm)[1] # get the name of the first taxa
  first_char <- strsplit(first_taxa_name,"")[[1]][1] # get the first character from the name of the first taxa
  if (first_char == "t"){
    # for phylo sims - sort cols and rows
    rn <- row.names(pdm) # get rownames
    rn <- substring(rn,2) #remove t from rownames
    rownames(pdm) <- sprintf("%03s",rn) # pad out rownames and reattach  (so will order numerically)
    cn <- colnames(pdm) # get colnames
    cn <- substring(cn,2)# remove t from colnames
    colnames(pdm) <- sprintf("%03s",cn) # pad out colnames and reattach (so will order numerically)
    pdm <- pdm[order(rownames(pdm)),order(colnames(pdm))] # order by number
  } else {
    # for SimBac sims - sort cols and rows
    rownames(pdm) <- sprintf("%03s",rownames(pdm)) # pad out rownames and reattach (so will order numerically)
    colnames(pdm) <- sprintf("%03s",colnames(pdm)) # pad out colnames and reattach (so will order numerically)
    pdm <- pdm[order(rownames(pdm)),order(colnames(pdm))] # order by number
  }
  # Then do upper triangle removal
  pdm[upper.tri(pdm)] <- 0 # Get upper triangle and replace upper triangle coordinates (TRUE) with 0
  return(pdm)
}

# Function to open a maximum likelihood distances file from IQ tree and return a pairwise distance matrix
mldist.pdm <- function(path){
  # Open the maximum likelihood distances file output from creating the tree in IQ-tree
  filename <- paste0(path,".mldist")
  pdm <- read.table(filename,header=FALSE,skip=1) # read in mldist file
  # After opening the mldist matrix, need to sort it so that the taxa are in the right order
  # So that when doing the subtractions, the pairwise distances will be coming from the same location (eg AB and AB, not BF and EH)
  # Place taxa in same order for both tree and alignment matrix
  # the first column will be the order of the taxa
  taxa <- as.character(pdm[,1])
  first_taxa_name <- taxa[1] # get the name of the first taxa
  first_char <- strsplit(first_taxa_name,"")[[1]][1] # get the first character from the name of the first taxa
  if (first_char == "t"){
    # Phylo sim - sort based on taxa naming conventions from the alignment
    pdm <- pdm[,2:ncol(pdm)] # remove col with names
    n <- substring(taxa,2) # remove the "t" character from the taxa name
    colnames(pdm) <- sprintf("%03s",n) # change the row names to the taxa labels, pad out colnames and reattach (so will order numerically)
    rownames(pdm) <- sprintf("%03s",n) # change the row names to the taxa labels, pad out rownames and reattach (so will order numerically)
    pdm <- pdm[order(rownames(pdm)),order(colnames(pdm))] # order by number
  } else {
    # otherwise is SimBac sim - sort based on taxa naming conventions from the alignment
    pdm <- pdm[,2:ncol(pdm)] # remove col with names
    colnames(pdm) <- sprintf("%03s",taxa) # change the row names to the taxa labels, pad out colnames and reattach (so will order numerically)
    rownames(pdm) <- sprintf("%03s",taxa) # change the row names to the taxa labels, pad out rownames and reattach (so will order numerically)
    pdm <- pdm[order(rownames(pdm)),order(colnames(pdm))] # order by number
  }
  # Then do upper triangle removal
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

# Get numerators and denominators for the test statistics you invented
get.TS.fractions <- function(path, iqpath, splitstree_path){
  # in case you forget: numerator goes on the top, denominator goes on the bottom
  # Get numbers from TS1
  tree_pdm <- iqtree.pdm(iqpath,path)
  ts1_num <- sum(tree_pdm)
  alignment_pdm <- mldist.pdm(path)
  ts1_denom <- sum(alignment_pdm)
  
  # Get fraction from TS4
  call.IQTREE(iqpath,path) # path = path to alignment
  # Calculate the split decomposition
  call.SplitsTree(splitstree_path,path,"split decomposition")
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
  # get the fraction details
  ts4_num <- tree_ii_sum
  ts4_denom <- all_ii_sum
  
  call.IQTREE(iqpath,path) # path = path to alignment
  # Calculate the split decomposition
  call.SplitsTree(splitstree_path,path,"neighbournet")
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
  ts5_num <- tree_ii_sum
  ts5_denom <- all_ii_sum
  
  # Collect all variables
  fracs <- c(ts1_num,ts1_denom,ts4_num,ts4_denom,ts5_num,ts5_denom)
  return(fracs)
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
  # average the difference matrix to get the test statistic
  # want the mean of the values in the upper triangle
  # extract all the values in the upper triangle of the diff_pdm matrix
  vals <- diff_pdm[lower.tri(diff_pdm)]
  # take the mean of all the values in the upper triangle to get the mean absolute difference pairwise distance
  ts <- mean(vals) 
  
  # return the test statistic (the mean calculated just above)
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

# Tree proportion
# Find which splits are in the tree and sum those split weights, divide by sum of all split weights#
tree.proportion <- function(iqpath, splitstree_path, path, network_algorithm = "neighbournet", trimmed = FALSE, tree_path = FALSE, run_IQTREE = FALSE){
  # network_algorithm takes values "split decomposition" or "neighbournet", default is "neighbournet"
  # trimmed takes values TRUE or FALSE, default is FALSE (all branches in network included in tree)
  
  # Only runs IQTree if run)IQTREE is set to TRUE
  # Then it will only run if it hasn't run before (won't run if .iqtree file exists)
  if (run_IQTREE == TRUE){
    # Run IQ-tree if it hasn't already been run
    call.IQTREE(iqpath,path) # path = path to alignment
  }
  
  # Calculate the split decomposition
  call.SplitsTree(splitstree_path,path,network_algorithm)
  # Retrieve the file name for the splits output file
  splits.filepath <- splits.filename(path)
  
  # Get the number of splits by opening the splits block and reading the number after "nsplits-"
  splits_text <- readLines(splits.filepath)
  splits_row <- splits_text[grep("nsplits",splits_text)]
  nsplits <- str_extract_all(splits_row,"[0-9]+")[[1]][2]
  # Check whether there are valid splits in the splits block
  # If there are not, return the tree proportion as NA
  # If there are, calculate the tree proportion and return it
  if (is.na(suppressWarnings(as.numeric(nsplits)))==TRUE){
    # If nsplits is not a number, return the ts as NA (can't calculate a tree proportion)
    ts <- NA
  } else if (as.numeric(nsplits) == 0){
    # If no splits, return NA (can't calculate a tree proportion)
    ts <- NA
  } else if (as.numeric(nsplits) > 0){
    # Otherwise, if the network does contain more than 0 splits:
    # Extract the splits 
    splits <- read.nexus.splits(splits.filepath) # Open the splits from SplitsTree
    # If no tree is provided, open the tree estimated by IQ_tree
    if (tree_path == FALSE){
      ## Open the tree estimated by IQ-TREE
      tree <- open.tree(path)
    } else {
      # if a tree is provided, open that tree
      tree <- read.tree(tree_path)
    }
    
    # Create empty sums to store the sum of isolation index weights in
    tree_ii_sum <- 0 
    all_ii_sum <- 0
    trimmed_ii_sum <- 0
    
    # Iterate through each of the rows in the splits dataframe and collect the splits
    for (i in 1:length(splits)){
      print(i)
      # For each split, get the split weight, whether it's trivial, and whether it's in the tree
      split_atts <- split.attributes(splits[i],tree)
      
      # Add the isolation index to the sum of all isolation indexes
      all_ii_sum <- all_ii_sum + split_atts$isolation_index
      # If the split is not trivial, add the isolation index to the sum of all non-trivial split isolation indexes
      if (split_atts$is_split_trivial == FALSE){
        trimmed_ii_sum <- trimmed_ii_sum + split_atts$isolation_index
      }
      
      # Add up the isolation index if the split is in the tree
      if (split_atts$is_split_in_tree == TRUE){
        if (trimmed == FALSE){
          # If not trimmed, add every monophyletic split up (include all trivial and non-trivial splits)
          tree_ii_sum <- tree_ii_sum + split_atts$isolation_index
        } else if (trimmed == TRUE){
          # If not trimmed, add ONLY the non-trivial splits (e.g. only splits with >1 species - ignore terminal branches)
          if (split_atts$is_split_trivial == FALSE){
            # Check that the split is non-trivial: if it is, add it to the running total
            tree_ii_sum <- tree_ii_sum + split_atts$isolation_index
          }
        }
      }
    }
    
    # Calculate the tree proportion by dividing the sum of split weights in the tree by the sum of all split weights
    # Use the tree_ii_sum which will already have pruned out the trivial splits if trimmed == TRUE
    if (trimmed == FALSE){
      # If no trimming, calculate by dividing by the sum of ALL isolation indexes
      ts <- tree_ii_sum / all_ii_sum
    } else if (trimmed == TRUE){
      # If trimming trivial splits, calculate test statistic by dividing by only the non-trivial splits
      ts <- tree_ii_sum / trimmed_ii_sum
    }
  }
  # Return the test statistic result
  return(ts)
}

# Function to get details about a split: whether it's monophyletic (and if so the isolation index), and whether it's trivial
split.attributes <- function(split, tree){
  # Extract the bipartition subsets from the tree
  ss1 <- split[[1]] # get the indices of all taxa in the split
  taxa <- attr(split,"labels") # get the names of all taxa in the tree
  ss1_taxa <- taxa[ss1] # use the split indices to get the taxa names from the split
  ss2_taxa <- setdiff(taxa,ss1_taxa) # take the elements from taxa not in subset 1 to get the elements in ss2
  # Test whether the split is trivial (if it is trivial, either ss1 or ss1 will be 1)
  if ((length(ss1_taxa) == 1) || (length(ss2_taxa) == 1)){
    trivial <- TRUE
  } else {
    trivial <- FALSE
  }
  # Remove any asterisks in the split labels to match the taxa names in the tree
  asterisk_taxa <- grep("\\*",taxa)
  if (length(asterisk_taxa)>0){
    # If any asterisk taxa exist, rename them by replacing the * with a _ - this is what happens when you open the IQTree tree
    taxa <- gsub("\\*","_",taxa)
    ss1_taxa <- taxa[ss1]
    ss2_taxa <- setdiff(taxa,ss1_taxa)
  }
  
  # Test whether the split is in the tree by seeing whether ss1 and ss2 are monophyletic
  ss1_mono <- is.monophyletic(tree,ss1_taxa)
  ss2_mono <- is.monophyletic(tree,ss2_taxa)
  # If both are true then add the split weight to the sum of split weights from splits in the tree
  # (if one is true, the other must be true - reduces to (A,B) where A and B are the two monophyletic clades)
  # Checks both just to be sure 
  # If they're not both tre
  if (ss1_mono == TRUE && ss2_mono == TRUE){
    in.tree <- TRUE # Split is in tree
  } else {
    in.tree <- FALSE # Split is not in tree
  }
  ii <- attr(split, "weights") # isolation index will equal the weight of the split (comes directly from the network) - this doesn't depend on whether the split is in the tree or not
  vals <- c(as.list(in.tree),ii,trivial)
  names(vals) <- c("is_split_in_tree","isolation_index","is_split_trivial")
  return(vals)
}

# Function to sum all weights from a splits nexus file
sum.all.ii <- function(splits){
  total_ii <- sum(attr(splits,"weights"))
  return(total_ii)
}

# Run split decomposition using SplitsTree
call.SplitsTree <- function(splitstree_path,alignment_path,network_algorithm){
  suffix <- tail(strsplit(alignment_path,"\\.")[[1]],1) # get the file extension from the filename
  if (suffix == "fasta"){
    # If the file is a fasta file, convert it to nexus file format
    data <- read.fasta(alignment_path) # read in the fasta data
    alignment_path_converted <- paste0(alignment_path,".nexus") # create a name for the nexus alignment (just the fasta alignment with a .nexus tacked on)
    write.nexus.data(data, file = alignment_path_converted,format = "dna",interleaved = TRUE, datablock = FALSE) # write the output as a nexus file)
    # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
    nexus <- readLines(alignment_path_converted) # open the new nexus file
    ind <- grep("BEGIN CHARACTERS",nexus)+2 # find which line
    nexus[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
    writeLines(nexus,alignment_path_converted) # output the edited nexus file
    plot_path <- gsub(".fasta.nexus","",alignment_path_converted)
  } else if (suffix == "nexus"){
    alignment_path_converted <- alignment_path
    plot_path <- gsub(".nexus","",alignment_path_converted)
  }
  else if (suffix == "nex"){
    alignment_path_converted <- alignment_path
    plot_path <- gsub(".nex","",alignment_path_converted)
  }
  output_path <- splits.filename(alignment_path) # create an output path
  if (network_algorithm == "split decomposition"){
    # Create the splitstree command
    splitstree_command <- paste0(splitstree_path, " -g -x 'OPEN FILE=", alignment_path_converted,"; ASSUME chartransform =Uncorrected_P HandleAmbiguousStates=Ignore Normalize=true; ASSUME disttransform=SplitDecomposition; SAVE FILE=", output_path," REPLACE=yes; QUIT'")
    # to export images of the networks (not working in soma)
    # splitstree_command <- paste0(splitstree_path, " -g -x 'OPEN FILE=", alignment_path_converted,"; ASSUME chartransform =Uncorrected_P HandleAmbiguousStates=Ignore Normalize=true; ASSUME disttransform=SplitDecomposition; SAVE FILE=", output_path," REPLACE=yes; UPDATE; EXPORTGRAPHICS format=SVG file=",paste0(plot_path,"_splitDecomposition.svg")," REPLACE=yes; UPDATE; EXPORTGRAPHICS format=PNG file=",paste0(plot_path,"_splitDecomposition.png")," REPLACE=yes; QUIT'")
  } else if (network_algorithm == "neighbournet"){
    splitstree_command <- paste0(splitstree_path, " -g -x 'OPEN FILE=", alignment_path_converted,"; ASSUME chartransform =Uncorrected_P HandleAmbiguousStates=Ignore Normalize=true; ASSUME disttransform=NeighbourNet; SAVE FILE=", output_path," REPLACE=yes; QUIT'")
    # to export images of the networks (not working in soma)
    # splitstree_command <- paste0(splitstree_path, " -g -x 'OPEN FILE=", alignment_path_converted,"; ASSUME chartransform =Uncorrected_P HandleAmbiguousStates=Ignore Normalize=true; ASSUME disttransform=NeighbourNet; SAVE FILE=", output_path," REPLACE=yes; UPDATE; EXPORTGRAPHICS format=SVG file=",paste0(plot_path,"_NeighbourNet.svg")," REPLACE=yes; UPDATE; EXPORTGRAPHICS format=png file=",paste0(plot_path,"_NeighbourNet.png")," REPLACE=yes; QUIT'")
  }
  # Call splitstree and do the split decomposition, save the results (overwrite any existing results)
  system(splitstree_command) # Call the splitstree command using the system 
}

# Function to create a file name for the SplitsTree output
splits.filename <- function(alignment_path){
  suffix <- tail(strsplit(alignment_path,"\\.")[[1]],1) # get the file extension from the filename
  # Create a identifiable name for the output file from splitstree
  if (suffix == "fasta"){
    # Add "_splits" into the filename (so can find file later)
    output_path <- gsub(".fasta","_splits.nexus",alignment_path)
  }
  if (suffix == "nexus"){
    # Add "_splits" into the filename (so can find file later)
    output_path <- gsub(".nexus","_splits.nexus",alignment_path)
  }
  if (suffix == "nex"){
    # Add "_splits" into the filename (so can find file later)
    output_path <- gsub(".nex","_splits.nex",alignment_path)
  }
  return(output_path)
}
