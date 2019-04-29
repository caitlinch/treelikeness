# Script to apply test statistics and parametric bootstrap to empirical data sets in the BenchmarkAlignments database

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

# Extract the file names of the alignments
files <- list.files(dirs,recursive = TRUE) # list all files
als <- paste0(dirs,files[grep(".nex",files)]) # get all the nexus files
als <- als[!als %in% als[grep(".nex.",als)]] # remove all non alignment files to leave only alignments

for (al in als){
  call.IQTREE(iqtree_path = "/Users/caitlincherryh/Documents/Honours/Executables/iqtree" , alignment_path = al)
}
