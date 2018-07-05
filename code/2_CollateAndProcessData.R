# R code to import and collate test statistic results, and to process the results

# Set working directory
maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # for work computer
setwd(maindir)

# Source files for functions
source(paste0(maindir,"code/func_process_data.R"))

# Collate data
external_df <- collate.csv("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments2/","external")
external_df <- simplify.SimBac(external_df)

internal_df <- collate.csv("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments2/","internal")
internal_df <- simplify.SimBac(internal_df)

phylo_df <- collate.csv("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments2/","pattern")
phylo_df <- simplify.phylo(phylo_df)