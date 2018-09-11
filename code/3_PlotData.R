# R code to plot simulation results
# Specify which file paths to use

run_location = "mac"
# run_location = "soma"

if (run_location == "mac"){
  # Set file paths etc
  input_folder <- "/Users/caitlincherryh/Documents/Results/collatedOutput/"
  # Set working directory
  maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # for work computer
} else if (run_location == "soma") {
  # Set file paths etc
  raw_data_folder <- "/data/caitlin/treelikeness/output/"
  output_folder <- "/data/caitlin/treelikeness/results/"
  
  # Set working directory
  maindir <- "/data/caitlin/treelikeness/"
}

# Source files for functions
source(paste0(maindir,"code/func_process_data.R"))

# load required libraries
library(ggplot2)
library(reshape2)
library(gridExtra)