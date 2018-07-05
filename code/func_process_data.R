# Functions to process and plot test statistics resulting from simulations

directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments2/"
id <- "CheckTestStats"
id <- "external"

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
  print(csv_paths)
  # Set the number of rows in the dataframe (will equal the number of folders - one folder per simulation)
  num_rows <- length(csv_paths)
  # Open all the csv files, store the results as a list
  output_list <- lapply(csv_paths,read.csv)
  
  # Create the empty dataframe
  output_df <- data.frame(matrix(nrow = num_rows, ncol = 32))
  lapply(1:num_rows, extract.one.csv, paths = csv_paths, df = output_df)
  return(output_df)
}

