# R functions to calculate the split decomposition for a pairwise distance matrix

library(TreeSim)
library(phytools)
library(seqinr)
library(ape)
library(phangorn)

# Get all options for n choose 2 from a given set
choose2 <- function(collection){
  if (length(collection) == 1){
    combinations <- list(c(collection[1],collection[1])) # If only one element, just send it back twice
  } else {
    # If more than 1 element, get all options for n choose 2
    combinations1 <- combn(collection, 2,simplify = FALSE) 
    # Also get all options for choosing each element twice
    combinations2 <- rep(list(),length(collection))
    for (i in 1:length(collection)){
      temp <- list(c(collection[i],collection[i]))
      combinations2[i] <- temp
    }
    # Need to combine the combination list and the choose-each-element-twice list
    combinations <- c(combinations1,combinations2)
  }
  return(combinations)
}


# Locate a pairwise distance in a matrix
# could change this to be matrix[a,b] and matrix[b,a] but user would have to submit a matrix with row and column names included
pd <- function(a,b,names,matrix){
  a_ind <- which(names == a) # get the indices
  b_ind <- which(names == b)
  a_pd <- matrix[a_ind,b_ind]  # Check to see whether the pairwise distance is stored under
  b_pd <- matrix[b_ind,a_ind ] # d[A,B] or d[B,A]
  dist <- max(a_pd,b_pd) # If both entries are 0, will have a pairwise distance of 0 - otherwise will be largest value
  return(dist)
}


# Calculate the inequality (whether the four-point condition holds) for a single quartet 
inequality <- function(i,j,k,l,names,matrix){
  d1 <- pd(i,j,names,matrix) + pd(k,l,names,matrix)
  d2 <- pd(i,k,names,matrix) + pd(j,l,names,matrix)
  d3 <- pd(i,l,names,matrix) + pd(j,k,names,matrix)
  if (d1 < (max(d2,d3))){
    # result = TRUE
    result <- 1 # If all quartets hold, will be 100% (a d-split)
  } else if (d1 >= max(d2,d3)){
    # result = FALSE
    result <- 0 # the more quartets don't hold the inequality, the lower the percentage will be 
  }
  return(result)
}


# Check whether one subset is a d-split or not
issplit <- function(partition,names,matrix,threshold){
  # Create a list to store 1 if the inequality holds and 0 if it doesn't
  # Used to check the percentage against the threshold (both must be 1 for a d-split as in Bandelt and Dress 1992)
  its <- c()
  # Find all combinations for the first subset
  combs1 <- choose2(partition[[1]])
  # Find all combinations for the second subset
  combs2 <- choose2(partition[[2]])
  # Iterate through first subset combinations
  for (x in combs1){
    # Iterate through second subset combinations
    for (y in combs2){
      # Calculate inequality for those 4 taxa
      i <- x[1]
      j <- x[2]
      k <- y[1]
      l <- y[2]
      calculate_inequality <- inequality(i,j,k,l,names,matrix)
      its <- c(its,calculate_inequality)
    }
  }
  summed <- sum(its) # sum all the iterations
  percentage <- summed/length(its) # divide the sum by the number of iterations to get a decimal value
  if (percentage >= threshold){
    # If inequality holds, sum will be equal to or greater than the threshold
    # For threshold = 1 (a pure d-split), the result will only be a split if all quartets hold the inequality
    # This is a split!
    split_result = TRUE
  } 
  if (percentage < threshold){
    # If doesn't hold, the percentage is lower than the threshold (e.g. some of te inequalities are false for a pure d-split)
    # This is not a split
    split_result = FALSE
  }
  return(split_result)
}


# Separate into partitions of two disjoint subsets
partition2 <- function(collection){
  # Need to go through the number of elements in each subset
  # n choose 1
  # n choose 2
  # ....
  # up to n choose (n/2) - for 6 need n choose 3, for 7 need n choose 3
  n <- length(collection)
  k <- (1:(floor(n/2)))
  combs <- c()
  for (i in k){
    # Get all options of n
    comb <- combn(collection, i, simplify = FALSE)
    combs <- c(combs,comb)
  }
  parts <- rep(list(),length(combs)) # initialise collection of partitions
  ind = 1 
  for (i in combs){
    indices_in <- which(collection %in% i)
    indices_out <- which(collection %in% setdiff(collection,i))
    part <- list(a = c(collection[indices_in]), b = c(collection[indices_out]))
    parts[[ind]] <- part 
    ind = ind+1
  }
  return(parts)
}


# Calculate the isolation index
isolation_index <- function(partition,names,matrix){
  # Isolation index formula (from Bandelt and Dress):
  # a = 1/2 min(max{d_ij+d_kl,d_ik+d_jl,d_il+d_jk}-d_ij-d_kl)
  mins <- c()
  # Find all combinations for the first subset
  combs1 <- choose2(partition[[1]])
  # Find all combinations for the second subset
  combs2 <- choose2(partition[[2]])
  # Try all quartets
  for (x in combs1){
    for (y in combs2){
      i = x[[1]]
      j = x[[2]]
      k = y[[1]]
      l = y[[2]]
      # Calculate all distances
      d1 <- pd(i,j,names,matrix) + pd(k,l,names,matrix)
      d2 <- pd(i,k,names,matrix) + pd(j,l,names,matrix)
      d3 <- pd(i,l,names,matrix) + pd(j,k,names,matrix)
      temp = max(d1,d2,d3) - d1
      mins = c(mins,temp)
    }
  }
  minimum = min(mins)
  alpha = 0.5*minimum
  return(alpha)
}


# Create a matrix of the right size
make.matrix <- function(num){
  mat <- matrix(0,nrow = num,ncol = num)
  return(mat)
}


# Create a matrix based on the partition
make.splitmatrix <- function(partition,names,alpha){
  # For a tree AB-CD, the distance AB and CD should be 0, and AC=BC=AD=BC=alpha (alpha is the isolation index)
  # The matrix would then be [[0 0 a a], [ 0,0,a,a], [a,a,0,0], [a,a,0,0]]
  # Create a matrix of the same number of rows and columns as the original matrix
  mat <- make.matrix(length(names))
  pd <- choose2(names)
  ss1 <- partition[[1]]
  ss2 <- partition[[2]]
  # For each possible combination of 2 taxa:
  for (i in pd){
    # If they're not in the same subset, change that distance to alpha (the isolation index)
    if ((i[1] %in% ss1 & i[2] %in% ss2) | (i[1] %in% ss2 & i[2] %in% ss1)){
      mat[which(names==i[2]),which(names==i[1])] <- alpha
    } 
  }
  return(mat)
}

# Calculate all partitions
# Specify threshold as a decimal (e.g. 1, 0.7, 0.5).
# If threshold = 1, all quartets must meet the four point condition for the split to be a d-split
# If threshold < 1, at least that percentage of quartets must meet the four point condition for the split to be counted as a split
# The threshold allows for a relaxed form of split decomposition 
split.decomposition <- function(taxa_names,distance_matrix,threshold = 1){
  sets = partition2(taxa_names)
  summed_matrix = make.matrix(length(taxa_names)) # make matrix using length of names vector
  initialise = TRUE
  for (part in sets){
    # Go through each partition and check if it's a split
    split = issplit(part,taxa_names,distance_matrix,threshold)
    if (split == TRUE){
      # If it is a split calculate the isolation index
      ii <- isolation_index(part,taxa_names,distance_matrix)
      # Use that isolation index to create a matrix
      matrix <- make.splitmatrix(part,taxa_names,ii)
      summed_matrix <- summed_matrix + matrix
      split_info <- list("partition" = part, "isolation_index" = ii, "matrix" = matrix) # Create a little list of information about the split
      # Add to the list of splits
      if (initialise == TRUE){
        splits_list <- c(list(split_info)) # initialise splits list if it doesn't exist
        initialise <- FALSE
      } else if (initialise == FALSE){
        splits_list <- c(splits_list,list(split_info)) # add to splits list if it does exist
      }
    }
  }
  # Create the output for the split decomposition
  op <- list(splits_list, summed_matrix)
  return(op)
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
  return(output_path)
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
  } else if (suffix == "nexus"){
    alignment_path_converted <- alignment_path
  }
  output_path <- splits.filename(alignment_path) # create an output path
  if (network_algorithm == "split decomposition"){
    # Create the splitstree command
    splitstree_command <- paste0(splitstree_path, " -g -x 'OPEN FILE=", alignment_path_converted,"; ASSUME chartransform =Uncorrected_P HandleAmbiguousStates=Ignore Normalize=true; ASSUME disttransform=SplitDecomposition; SAVE FILE=", output_path," REPLACE=yes; QUIT'")
  } else if (network_algorithm == "neighbournet"){
    splitstree_command <- paste0(splitstree_path, " -g -x 'OPEN FILE=", alignment_path_converted,"; ASSUME chartransform =Uncorrected_P HandleAmbiguousStates=Ignore Normalize=true; ASSUME disttransform=NeighbourNet; SAVE FILE=", output_path," REPLACE=yes; QUIT'")
  }
  # Call splitstree and do the split decomposition, save the results (overwrite any existing results)
  system(splitstree_command) # Call the splitstree command using the system 
}

