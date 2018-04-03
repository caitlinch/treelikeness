# Code to run the parametric bootstraps for the relevant datasets

# Open packages
library(ape)
library(TreeSim)
library(phytools)
library(phangorn)
library(base)

# Set working directory
# maindir <- "/Users/caitlin/Repositories/treelikeness/" # for laptop
maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # for work computer
setwd(maindir)

# Source files for functions
source(paste0(maindir,"code/split_decomposition.R"))
source(paste0(maindir,"code/parametric_bootstrap.R"))
source(paste0(maindir,"code/test_statistic.R"))

##### Code from first drafts of functions ####
## Test code for split decomposition functions
# # Enter distance matrix and taxa labels for practicing
# d <- t(matrix(c(0,0,0,0,0,0,0,4.0,0,0,0,0,0,0,5.0,1.0,0,0,0,0,0,7.0,3.0,2.0,0,0,0,0,13.0,9.0,8.0,6.0,0,0,
#                 0,8.0,12.0,13.0,11.0,5.0,0,0,6.0,10.0,11.0,13.0,7.0,2.0,0),
#               nrow = 7, ncol = 7)) # transpose as R filles by columns first not by rows first
# taxa <- c("A","B","C","D","E","F","G")
# rownames(d) <- taxa # label rownames with taxa
# colnames(d) <- taxa # label colnames with taxa
# a = split_decomposition(taxa,d)


# Code from writing test statistics to call IQ-TREE on each alignment  
# ## Test code to call test statistics and run them
# # # Input variables and files
# alignment_path <- "/Users/caitlin/Repositories/treelikeness/raw_data" # folder where alignment is located
# alignment_paths <- list.dirs(alignment_path)
# alignment_paths <- paste0(alignment_paths[2:length(alignment_paths)],"/") # to run all alignments in directory
# alignment_file <- "alignment.nex" # name of alignment 
# 
# for (alignment in alignment_paths){
#    print(alignment)
#    system(paste0(iqtree_path," -s ",alignment,alignment_file," -nt AUTO -redo"))
# }


##### Bootstrap code ####

### User Input Parameters
nbootstrap        <- 3 # the number of parametric bootstraps you'd like to run
alignments_folder <- paste0(maindir,"raw_data/") # The folder containing all of the datasets of interest
iqtree_path       <- "/Applications/iqtree/bin/iqtree" # location of IQ-tree program 


### To run the parametric bootstrap nbootstrap (e.g. 9/99/999) times on each alignment in the alignments folder:
# get a list of the folders (each folder = 1 alignment) in the raw_data folder
alignments <- paste0(alignments_folder,list.files(alignments_folder))
# iterate through each of the alignments
for (al in alignments){
  # extract information from the alignment.nex summary text file and the .iqtree file in the al folder
  # need these files to extract information to run the parametric bootstrap
  # parameters in order: number of taxa, number of sites, substitution database
  params <- get.simulation.parameters(al)
  # generate a list of ids for each simulation
  ids <- sprintf("%04d",1:nbootstrap)
  # create somewhere to store all the test statistics from the bootstraps
  all_ts <- c()
  for (id in ids){
    bootstrap_folder <- paste0(al,"/",id,"/") # make the filepath this bootstrap will be stored in
    # create a folder to store all the information about this bootstrap (if it doesn't already exist)
    if (dir.exists(bootstrap_folder) == FALSE){
      dir.create(bootstrap_folder)
    }
    # call the parametric bootstrap function on the folder of interest
    # pick test statistic based on numbering from grant proposal (Cherryh 2018)
    # test statistic = 1 calls the divide matrix test statistic
    # test statistic = 2 calls the matrix differences test statistic (unfinished)
    # test statistic = 3 calls the split decomposition test statistic (unfinished)
    bs_ts <- do.1.bootstrap(iqtree_path,bootstrap_folder,params,test_statistic=1)
    all_ts <- c(all_ts,bs_ts)
  }
  # For the alignment, call the test statistic
  ts <- c() # Put function call here
  # Save the actual test statistic
  all_ts <- c(ts,all_ts)
  # Create a list of names for ids for the output table
  names <- c("alignment",sprintf("bootstrap_%04d",1:(length(all_ts)-1))) # create a vector of names
  # Create an output data frame and label the columns
  ts_df <- data.frame(names,all_ts)
  names(ts_df) <- c("id","test_statistic")
  # Save the output dataframe
  ts_df_filename <- paste0(al,"/bootstrap_test_statistics_",as.character(Sys.Date),".csv")
  write.csv(ts_df,name = ts_df_filename)
}









