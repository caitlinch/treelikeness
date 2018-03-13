# R program to calculate a test statistic for quantifying tree likeness

# Input variables and files
alignment_path <- "/Users/caitlin/Repositories/treelikeness/raw_data/Crawford_2012/" # folder where alignment is located
alignment_file <- "alignment.nex" # name of alignment 
iqtree_path <- "/Applications/iqtree/bin/iqtree" # location of IQ-tree program

# Open packages
library(ape)

test_statistic <- function(iqpath,path,file){
    # Given an alignment, get a tree from IQ-tree and find the sum of the pairwise distance matrix
    # system(paste0(iqpath," -s ",alignment_path,alignment_file," -nt AUTO -redo")) # call IQ-tree!
    tree <- read.tree(paste0(path,file,".treefile")) # read in the tree created in IQ-tree
    tree_pdm <- cophenetic.phylo(tree) # get pairwise distance matrix for the taxa - note that this mirrors across the diagonal! (doesn't have half filled with 0's)
    tree_pdm[upper.tri(tree_pdm)] <- 0 # Get upper triangle and replace upper triangle coordinates (TRUE) with 0
    print(tree_pdm) # to check
    tree_sum <- sum(tree_pdm) # sum up all the pairwise distances
    tree_pdm <- tree_pdm/tree_sum # scale matrix entries to add to 1
    print(tree_pdm) # print scaled pdm matrix
    tree_sum <- sum(tree_pdm) # get new tree sum
    print(tree_sum) # print tree sum
    
    # Calculate pairwise distances from the alignment
    # Open the maximum likelihood distances file output when creating the tree in IQ-tree
    mldist_file <- paste0(path,file,".mldist")
    mldist_pdm <- read.table(mldist_file,header=FALSE,skip=1) # read in mldist file
    mldist_pdm <- mldist_pdm[2:ncol(mldist_pdm)] # remove first column (taxa names)
    mldist_pdm[upper.tri(mldist_pdm)] <- 0 # replace upper triangle coordinates (TRUE) with 0
    print(mldist_pdm) # to check
    mldist_sum <- sum(mldist_pdm) # sum all the pairwise distances
    mldist_pdm <- mldist_pdm/mldist_sum # scale matrix entries to add to 1
    print(mldist_pdm) # print scaled pdm matrix
    mldist_sum <- sum(mldist_pdm) # get new tree sum
    print(mldist_sum) # print tree sum
    
    # Divide sums
    ts <- tree_sum/mldist_sum
    divisor_pdm <- (tree_pdm/mldist_pdm)
    print(divisor_pdm)
    print(sum(divisor_pdm,na.rm = TRUE))
    return()
}