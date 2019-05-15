# Script to apply test statistics and parametric bootstrap to empirical data sets in the BenchmarkAlignments database

library(ape)
library(parallel)
library(phangorn)
library(phytools)
library(seqinr)
library(stringr)
library(TreeSim)

# run_location = "mac"
run_location = "soma"

if (run_location == "mac"){
  BA_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/BA_testSet/"
  output_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/BA_testSet_output/"
  maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # where the code is
  exec_folder <- "/Users/caitlincherryh/Documents/Honours/Executables/"
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("3seq","iqtree","Phi","SimBac","SplitsTree.app/Contents/MacOS/JavaApplicationStub")
  exec_paths <- paste0(exec_folder,exec_paths)
  names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")
  source(paste0(maindir,"code/func_BA.R"))
} else if (run_location=="soma"){
  BA_dir <- "/data/caitlin/treelikeness/BenchmarkAlignments_DataSubSet/"
  output_dir <- "/data/caitlin/treelikeness/BenchmarkAlignments_DataSubSet_Results/"
  maindir <- "/data/caitlin/treelikeness/" # where the code is
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/data/caitlin/linux_executables/3seq/3seq","/data/caitlin/linux_executables/iqtree/bin/iqtree","/data/caitlin/linux_executables/PhiPack/Phi",
                  "/data/caitlin/linux_executables/SimBac/SimBac","/data/caitlin/splitstree4/SplitsTree")
  names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")
  source(paste0(maindir,"code/func_BA_parallel.R")) # run code parallel
}

# Source files for functions
source(paste0(maindir,"code/func_test_statistic.R"))
source(paste0(maindir,"code/func_process_data.R"))
source(paste0(maindir,"code/func_parametric_bootstrap.R"))

# Extract the file names of the alignments
files <- list.files(BA_dir,recursive = TRUE) # list all files
als <- paste0(BA_dir,files[grep(".nex",files)]) # get all the nexus files
als <- als[!als %in% als[grep(".nex.",als)]] # remove all non alignment files to leave only alignments
als <- als[!als %in% als[grep("bootstrapReplicate",als)]] # remove all bootstrap alignments (if any present)
als <- als[!als %in% als[grep("alignment.nex",als)]] # remove full alignments (only want to run per loci)
als <- sort(als,decreasing = TRUE) # reverse list so it runs the HUGE mammal dataset first

# Calculate the test statistics and run the bootstraps
# To run for one alignment: empirical.bootstraps.wrapper(empirical_alignment_path = empirical_alignment_path, program_paths = program_paths, number_of_replicates = 9)
if (run_location=="soma"){
  lapply(als,empirical.bootstraps.wrapper, program_paths = exec_paths, number_of_replicates = 199) 
} else if (run_location=="mac"){
  lapply(als,empirical.bootstraps.wrapper, program_paths = exec_paths, number_of_replicates = 9) 
}

# Collate all the results
results_file <- paste0(output_dir,basename(BA_dir),"_completeResults.csv")
df <- collate.bootstraps(directory = BA_dir, file.name = "pValues", id = "", output.file.name = results_file)

