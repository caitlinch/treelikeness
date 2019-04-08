# Source functions
maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/"
source(paste0(maindir,"code/func_split_decomposition.R"))
source(paste0(maindir,"code/func_test_statistic.R"))
source(paste0(maindir,"code/func_create_alignments.R"))
source(paste0(maindir,"code/func_process_data.R"))
source(paste0(maindir,"code/func_parametric_bootstrap.R"))
tree_folder <- paste0(maindir,"trees/")

# Set filepaths
al_path <- "/Users/caitlincherryh/Documents/Honours/TestAlignmentResults/5_newTS/divergent_alignment.nexus"
results_folder <- "/Users/caitlincherryh/Documents/Honours/TestAlignmentResults/5_newTS/"
main_dir <- "/Users/caitlincherryh/Documents/Honours/TestAlignmentResults/5_newTS/"
exec_folder <- "/Users/caitlincherryh/Documents/Honours/Executables/"
exec_paths <- c("3seq","iqtree","Phi","SimBac","SplitsTree.app/Contents/MacOS/JavaApplicationStub")
exec_paths <- paste0(exec_folder,exec_paths)
names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")

# Call decomposition statistic
x <- tree.proportion(iqpath = exec_paths[["IQTree"]], splitstree_path = exec_paths[["SplitsTree"]], path = al_path, network_algorithm = "neighbournet", trimmed = FALSE)
y <- tree.proportion(iqpath = exec_paths[["IQTree"]], splitstree_path = exec_paths[["SplitsTree"]], path = al_path, network_algorithm = "neighbournet", trimmed = TRUE)