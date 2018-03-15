# R functions to calculate the split decomposition for a pairwise distance matrix

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
    result <- 0
  } else if (d1 >= max(d2,d3)){
    # result = FALSE
    result <- 1
  }
  return(result)
}


# Check whether one subset is a d-split or not
issplit <- function(partition,names,matrix){
  # Create a list to store 0 if the inequality holds and 1 if it doesn't
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
  summed <- sum(its)
  if (summed == 0){
    # If inequality holds, sum will be 0
    split_result = TRUE
  } 
  if (summed > 0){
    # If doesn't hold, sum larger than 0
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
make_matrix <- function(num){
  mat <- matrix(0,nrow = num,ncol = num)
  return(mat)
}


# Create a matrix based on the partition
make_splitmatrix <- function(partition,names,alpha){
  # For a tree AB-CD, the distance AB and CD should be 0, and AC=BC=AD=BC=alpha (alpha is the isolation index)
  # The matrix would then be [[0 0 a a], [ 0,0,a,a], [a,a,0,0], [a,a,0,0]]
  # Create a matrix of the same number of rows and columns as the original matrix
  mat <- make_matrix(length(names))
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
split_decomposition <- function(taxa_names,distance_matrix){
  sets = partition2(taxa_names)
  summed_matrix = make_matrix(length(taxa_names)) # make matrix using length of names vector
  splits_list <- c() # initialise collection of splits
  for (part in sets){
    # Go through each partition and check if it's a split
    split = issplit(part,taxa_names,distance_matrix)
    if (split == TRUE){
      # If it is a split calculate the isolation index
      ii <- isolation_index(part,taxa_names,distance_matrix)
      # Use that isolation index to create a matrix
      matrix <- make_splitmatrix(part,taxa_names,ii)
      summed_matrix <- summed_matrix + matrix
      split_info <- list("partition" = part, "isolation index" = ii, "matrix" = matrix) # Create a little list of information about the split
      print(split_info)
      # Add to the list of splits
      splits_list <- append(splits_list,split_info)
    }
  }
  # Create the output for the split decomposition
  op <- list(splits_list, summed_matrix)
  return(op)
}



# Enter distance matrix and taxa labels for practicing
d <- t(matrix(c(0,0,0,0,0,0,0,4.0,0,0,0,0,0,0,5.0,1.0,0,0,0,0,0,7.0,3.0,2.0,0,0,0,0,13.0,9.0,8.0,6.0,0,0,
                0,8.0,12.0,13.0,11.0,5.0,0,0,6.0,10.0,11.0,13.0,7.0,2.0,0),
              nrow = 7, ncol = 7)) # transpose as R filles by columns first not by rows first
taxa <- c("A","B","C","D","E","F","G")
rownames(d) <- taxa # label rownames with taxa
colnames(d) <- taxa # label colnames with taxa
a = split_decomposition(taxa,d)

