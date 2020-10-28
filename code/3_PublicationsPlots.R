# R code to import dataframes of test statistic and statistical test results and create plots to summarise the results
# Sourcing this file will open five dataframes and output a number of plots
# Final result is a number of plots displaying test statistic values under perturbation of various simulation factors

##### Step 1: Open packages #####
library(ggplot2)
library(ggpmisc)
library(patchwork)




##### Step 2: Uncomment and set the file paths for output folders, executables, and the identifying name for this run #####
# results_folder <- the folder where the result csvs will be placed (I use same results_folder in Parts 1 - 4.)
# plots_folder <- the folder where the plots will be stored
# maindir <- "treelikeness" repository location
# run_id <- if "run.id  = FALSE", program extracts run_id from input parameter file names 
#        <- otherwise, run_id will be set to whatever the user inputs here (e.g. "run_id = 'replicateAnalysis' ")

# results_folder <- ""
# plots_folder <- ""
# maindir <- ""
# run_id = FALSE
# tree_length = 0.5  # choose a value for total tree depth: 0.05, 0.10, 0.5, or 1
# this value will be used as tree depth in all plots EXCEPT the plots that contain all 4 tree depths for comparison


#__________________________________________Caitlin's paths (delete these if you're not Caitlin)______________________________________
results_folder <- "/Users/caitlincherryh/Documents/Honours/Results/simulations_20200304/output/"
plots_folder <- "/Users/caitlincherryh/Documents/Honours/Results/simulations_20200304/plots/"
maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/"
run_id = FALSE
tree_length = 0.5
#____________________________________________________________________________________________________________________________________



##### Step 3: Source function files #####
source(paste0(maindir,"code/func_process_data.R"))
# Extract run.id from the results folder name so the whole analysis has the same run.id
if (run_id == "FALSE"){
  run_id <- extract.run.id(results_folder)
}



##### Step 4: Open dataframes #####
# List all the files in the results_folder
melt_files <- list.files(results_folder)[grep("melted",list.files(results_folder))]
# Using the list of files in the results_folder, extract the melted csvs for each experiment (they are ready to plot!)
# p1:p3 are the dataframes for the test statistics calculated and estimated in experiments 1:3 respectively
ts1_df <- read.csv(paste0(results_folder,melt_files[grep("exp1_testStatistics",melt_files)]), stringsAsFactors = FALSE)
ts2_df <- read.csv(paste0(results_folder,melt_files[grep("exp2_testStatistics",melt_files)]), stringsAsFactors = FALSE)
ts3_df <- read.csv(paste0(results_folder,melt_files[grep("exp3_testStatistics",melt_files)]), stringsAsFactors = FALSE)
# bootstrap dataframe contains information about the p values (obtained for tree proportion using a parametric bootstrap)
# bs2:bs3 are the dataframes for the p values calculated and estimated in experiments 2:3 respectively
bs2_df <- read.csv(paste0(results_folder,melt_files[grep("exp2_p_value",melt_files)]), stringsAsFactors = FALSE)
bs3_df <- read.csv(paste0(results_folder,melt_files[grep("exp3_p_value",melt_files)]), stringsAsFactors = FALSE)


