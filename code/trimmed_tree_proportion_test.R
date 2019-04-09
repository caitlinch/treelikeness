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
al2_path <- "/Users/caitlincherryh/Documents/Honours/TestAlignmentResults/5_newTS/32_8event.fasta"
results_folder <- "/Users/caitlincherryh/Documents/Honours/TestAlignmentResults/5_newTS/"
main_dir <- "/Users/caitlincherryh/Documents/Honours/TestAlignmentResults/5_newTS/"
exec_folder <- "/Users/caitlincherryh/Documents/Honours/Executables/"
exec_paths <- c("3seq","iqtree","Phi","SimBac","SplitsTree.app/Contents/MacOS/JavaApplicationStub")
exec_paths <- paste0(exec_folder,exec_paths)
names(exec_paths) <- c("3seq","IQTree","Phi","SimBac","SplitsTree")

# Call decomposition statistic
x <- tree.proportion(iqpath = exec_paths[["IQTree"]], splitstree_path = exec_paths[["SplitsTree"]], path = al2_path, network_algorithm = "neighbournet", trimmed = FALSE)
y <- tree.proportion(iqpath = exec_paths[["IQTree"]], splitstree_path = exec_paths[["SplitsTree"]], path = al2_path, network_algorithm = "neighbournet", trimmed = TRUE)

# Call delta plot function
mldm <- mldist.pdm(al2_path)
pdmm <- as.matrix(mldm)
deltaplot_results <- delta.plot(pdmm,k = 101, plot = FALSE) # calculate the delta.plot
counts <- deltaplot_results$counts
intervals <- seq(0,1,(1/(length(counts)-1)))
deltaplot_df <- data.frame(intervals,counts)
names(deltaplot_df) <- c("intervals","counts")
deltaplot_df_name <- phylo.fixedtrees.output.folder(row)[4]
write.csv(df,file = deltaplot_df_name)
# Want to calculate the mean and median delta q value - unfortunately the delta.plot function doesn't output raw data, so make a pseudo data set using the histogram values
mean_dq <- mean(rep(deltaplot_df$intervals,deltaplot_df$counts)) # turn the interval data into a long list of "raw" values and calculate the mean
median_dq <- median(rep(deltaplot_df$intervals,deltaplot_df$counts)) # turn the interval data into a long list of "raw" values and calculate the median
