# Script to apply test statistics and parametric bootstrap to empirical data sets in the BenchmarkAlignments database

library(ape)
library(phangorn)

run_location = "mac"
# run_location = "soma"

if (run_location == "mac"){
  BA_dir <- "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/BA_testSet/"
  maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/"
  exec_folder <- "/Users/caitlincherryh/Documents/Honours/Executables/"
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("3seq","iqtree","Phi","SimBac","SplitsTree.app/Contents/MacOS/JavaApplicationStub")
  exec_paths <- paste0(exec_folder,exec_paths)
  names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")
} else if (run_location=="soma"){
  BA_dir <- ""
  maindir <- "/data/caitlin/treelikeness/"
  exec_paths <- c("/data/caitlin/linux_executables/3seq/3seq","/data/caitlin/linux_executables/iqtree/bin/iqtree","/data/caitlin/linux_executables/PhiPack/Phi",
                  "/data/caitlin/linux_executables/SimBac/SimBac","/data/caitlin/splitstree4/SplitsTree")
  names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")
}

# Source files for functions
source(paste0(maindir,"code/func_test_statistic.R"))
source(paste0(maindir,"code/func_parametric_bootstrap.R"))
source(paste0(maindir,"code/func_BA.R"))

# Extract the file names of the alignments
files <- list.files(BA_dir,recursive = TRUE) # list all files
als <- paste0(BA_dir,files[grep(".nex",files)]) # get all the nexus files
als <- als[!als %in% als[grep(".nex.",als)]] # remove all non alignment files to leave only alignments

alignment_path <- als[1]
program_paths <- exec_paths

tree.proportion(iqpath = program_paths[["IQTree"]], splitstree_path = program_paths[["SplitsTree"]],path = "/Users/caitlincherryh/Documents/Chapter01_TestStatistics_BenchmarkAlignments/BA_testSet/Anderson_2013/16S/16S_withTaxaBlock.nexus", network_algorithm = "split decomposition", trimmed = FALSE)

tree.proportion(iqpath = program_paths[["IQTree"]], splitstree_path = program_paths[["SplitsTree"]],path = "/Users/caitlincherryh/Documents/Honours/Executables/alignment.nexus", network_algorithm = "split decomposition", trimmed = FALSE)

#empirical.runTS(als[1],exec_paths)

#for (al in als[1]){
#  empirical.runTS(al,exec_paths)
#}
