# Quick script to extract the fractions for test statistics 1, 4 and 5 from the alignments to graph them

folder_paths <- list.files("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments3/")
# Get the positions at which the simulations containing the id are at
inds <- grep("internal", folder_paths)
# Get only the relevant folders (folders containing the id)
csv_paths <- folder_paths[inds]
# Make the path to each output csv file
csv_paths <- paste0("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments3/", csv_paths, "/alignment.fasta")
data <- lapply(csv_paths,get.TS.fractions,iqpath = iqpath, splitstree_path = splitstree_path)
fracs <- as.data.frame(matrix(unlist(data), ncol = 6, byrow = TRUE))
names(fracs)<- c("ts1_num","ts1_denom","ts4_num","ts4_denom","ts5_num","ts5_denom")
file_name <- paste0("/Users/caitlincherryh/Documents/TestAlignmentResults/0_prelim_mk3/","internal_testStatisticFractions_","fixedPDM",".csv")
write.csv(fracs, file = file_name, row.names = FALSE)

folder_paths <- list.files("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments3/")
# Get the positions at which the simulations containing the id are at
inds <- grep("external", folder_paths)
# Get only the relevant folders (folders containing the id)
csv_paths <- folder_paths[inds]
# Make the path to each output csv file
csv_paths <- paste0("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments3/", csv_paths, "/alignment.fasta")
data <- lapply(csv_paths,get.TS.fractions,iqpath = iqpath, splitstree_path = splitstree_path)
fracs <- as.data.frame(matrix(unlist(data), ncol = 6, byrow = TRUE))
names(fracs)<- c("ts1_num","ts1_denom","ts4_num","ts4_denom","ts5_num","ts5_denom")
file_name <- paste0("/Users/caitlincherryh/Documents/TestAlignmentResults/0_prelim_mk3/","external_testStatisticFractions_","fixedPDM",".csv")
write.csv(fracs, file = file_name, row.names = FALSE)

folder_paths <- list.files("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments3/")
# Get the positions at which the simulations containing the id are at
inds <- grep("2trees", folder_paths)
# Get only the relevant folders (folders containing the id)
csv_paths <- folder_paths[inds]
# Make the path to each output csv file
csv_paths <- paste0("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments3/", csv_paths, "/alignment.nexus")
data <- lapply(csv_paths,get.TS.fractions,iqpath = iqpath, splitstree_path = splitstree_path)
fracs <- as.data.frame(matrix(unlist(data), ncol = 6, byrow = TRUE))
names(fracs)<- c("ts1_num","ts1_denom","ts4_num","ts4_denom","ts5_num","ts5_denom")
file_name <- paste0("/Users/caitlincherryh/Documents/TestAlignmentResults/0_prelim_mk3/","2trees_testStatisticFractions_","_fixedPDM",".csv")
write.csv(fracs, file = file_name, row.names = FALSE)
