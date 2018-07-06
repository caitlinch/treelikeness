# Functions to process and plot test statistics resulting from simulations

# Function to collect all test statistic csv folders given an id to look for and a directory to look in
collate.csv <- function(directory,id){
  # Collect all the folders within the directory
  folder_paths <- list.files(directory)
  # Get the positions at which the simulations containing the id are at
  inds <- grep(id, folder_paths)
  # Get only the relevant folders (folders containing the id)
  csv_paths <- folder_paths[inds]
  # Make the path to each output csv file
  csv_paths <- paste0(directory, csv_paths, "/testStatistics.csv")
  # Set the number of rows in the dataframe (will equal the number of folders - one folder per simulation)
  num_rows <- length(csv_paths)
  # Open all the csv files, store the results as a list
  output_list <- lapply(csv_paths, read.csv, stringsAsFactors = FALSE)
  # Reduce the dataframe from a list into a matrix
  output_df <- Reduce(rbind, output_list)
  # Remove the "X" column (row number for smaller csv files)
  output_df <- output_df[,4:33]
  return(output_df)
}

# Extract the relevant columns for analysing the SimBac simulations
simplify.SimBac <- function(df){
  cols <- c(1:5,14:30) 
  df <- df[cols]
  return(df)
}

# Extract the relevant columns for analysing the phylogenetic simulations
simplify.phylo <- function(df){
  cols <- c(1:2,6:11,14:30) 
  df <- df[cols]
  return(df)
}

