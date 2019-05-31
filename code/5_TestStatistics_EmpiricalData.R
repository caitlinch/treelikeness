# Script to apply test statistics and parametric bootstrap to empirical data sets in the BenchmarkAlignments database

print("opening packages")
library(ape)
library(parallel)
library(phangorn)
library(phytools)
library(seqinr)
library(stringr)
library(TreeSim)

print("initialising namespace")

# run_location = "mac"
run_location = "soma"
# run_location = "nci"

if (run_location == "mac"){
  BA_dir <- "/Users/caitlincherryh/Documents/Repositories/BenchmarkAlignments_DataSubSet/"
  output_dir <- "/Users/caitlincherryh/Documents/Repositories/BenchmarkAlignments_DataSubSet/"
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
} else if (run_location == "nci"){
  BA_dir <- "/short/xf1/cac599/BenchmarkAlignments_DataSubSet/"
  output_dir <- "/short/xf1/cac599/BenchmarkAlignments_DataSubSet_Results/"
  maindir <- "/home/599/cac599/treelikeness/" # where the code is
  # Create a vector with all of the executable file paths
  # To access a path: exec_paths[["name"]]
  exec_paths <- c("/home/599/cac599/treelikeness/executables/3seq","/home/599/cac599/treelikeness/executables/iqtree","/home/599/cac599/treelikeness/executables/Phi",
                  "/home/599/cac599/treelikeness/executables/SimBac","/home/599/cac599/treelikeness/splitstree4/SplitsTree")
  names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")
  # source(paste0(maindir,"code/func_BA_parallel.R")) # run code parallel
  source(paste0(maindir,"code/func_BA.R"))
}

# Source files for functions
source(paste0(maindir,"code/func_test_statistic.R"))
source(paste0(maindir,"code/func_process_data.R"))
source(paste0(maindir,"code/func_parametric_bootstrap.R"))

# Extract the file names of the alignments
print("extracting alignments")
files <- list.files(BA_dir,recursive = TRUE) # list all files
als <- paste0(BA_dir,files[grep(".nex",files)]) # get all the nexus files
als <- als[!als %in% als[grep(".nex.",als)]] # remove all non alignment files to leave only alignments
als <- als[!als %in% als[grep("bootstrapReplicate",als)]] # remove all bootstrap alignments (if any present)
als <- als[!als %in% als[grep("alignment.nex",als)]] # remove full alignments (only want to run per loci)
# Nucleotide datasets ordered from smallest to largest number of partitions in original dataset
partition_order <- c("Worobey_2014h","Worobey_2014a","Worobey_2014f","Worobey_2014g","Worobey_2014e","Worobey_2014b",
                  "Worobey_2014c","Worobey_2014d","Looney_2016","Anderson_2013","Devitt_2013","Seago_2011","Siler_2013",
                  "Cognato_2001","Brown_2012","Wood_2012","Bergsten_2013","Murray_2013","Kawahara_2013","Sauquet_2011",
                  "Day_2013","Sharanowski_2011","Tolley_2013","Dornburg_2012","Unmack_2013","Rightmyer_2013","Horn_2014",
                  "Wainwright_2012","Near_2013","Pyron_2011","Oaks_2011","Lartillot_2012","Broughton_2013","Reddy_2017",
                  "Fong_2012","Cannon_2016_dna","Faircloth_2013","Moyle_2016","Leache_2015","Branstetter_2017",
                  "Crawford_2012","Smith_2014","Meiklejohn_2016","McCormack_2013","Richart_2015","Prebus_2017",
                  "Ran_2018_dna","Wu_2018_dna")
# Nucleotide datasets ordered from smallest to largest number of taxa
taxa_order <- c('Richart_2015','Smith_2014','Crawford_2012','Leache_2015','Meiklejohn_2016','Faircloth_2013',
                'McCormack_2013','Wood_2012','Bergsten_2013','Ran_2018_dna','Brown_2012','Cognato_2001',
                'Dornburg_2012','Prebus_2017','Sauquet_2011','Broughton_2013','Siler_2013','Devitt_2013',
                'Kawahara_2013','Cannon_2016_dna','Lartillot_2012','Oaks_2011','Wu_2018_dna','Worobey_2014c',
                'Rightmyer_2013','Seago_2011','Moyle_2016','Fong_2012','Unmack_2013','Sharanowski_2011',
                'Anderson_2013','Worobey_2014a','Day_2013','Branstetter_2017','Wainwright_2012','Horn_2014',
                'Tolley_2013','Reddy_2017','Murray_2013','Worobey_2014g','Worobey_2014b','Worobey_2014f',
                'Worobey_2014e','Worobey_2014h','Worobey_2014d','Near_2013','Looney_2016','Pyron_2011')
# Sort the list of alignments by the desired order
print("sorting alignments")
als_ordered <- c()
for (ds in taxa_order){
  als_ordered <- c(als_ordered,als[grep(ds,als)])
}
als <- als_ordered

# Calculate the test statistics and run the bootstraps
print("starting analysis")
# To run for one alignment: empirical.bootstraps.wrapper(empirical_alignment_path = empirical_alignment_path, program_paths = program_paths, number_of_replicates = 9)
if (run_location=="soma"){
  lapply(als,empirical.bootstraps.wrapper, program_paths = exec_paths, number_of_replicates = 99) 
} else if (run_location=="mac"){
  lapply(als,empirical.bootstraps.wrapper, program_paths = exec_paths, number_of_replicates = 9) 
}

# Collate all the results
results_file <- paste0(output_dir,basename(BA_dir),"_completeResults.csv")
df <- collate.bootstraps(directory = BA_dir, file.name = "pValues", id = "", output.file.name = results_file)

