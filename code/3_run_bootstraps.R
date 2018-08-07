# Code to run the parametric bootstraps for the relevant datasets

## User input 
source = "mac"

if (source=="mac"){
  test_alignment_folder <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testBootstrap1/Phylo_20_1300_NA_NA_NA_1_0.1_0.47_2trees_5/"
  iq_path <- "/Users/caitlincherryh/Documents/Executables/iqtree"
  splitstree_path <- "/Users/caitlincherryh/Documents/Executables/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
  n_reps <- 2 
}

## Main Code 
# Open packages
library(ape)
library(TreeSim)
library(phytools)
library(phangorn)
library(base)

# Source files for functions
source(paste0(maindir,"code/func_split_decomposition.R"))
source(paste0(maindir,"code/func_parametric_bootstrap.R"))
source(paste0(maindir,"code/func_test_statistic.R"))



