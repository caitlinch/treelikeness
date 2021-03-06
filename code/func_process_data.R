# Functions to process and plot test statistics resulting from simulations

# Function to collect all test statistic csv folders given an id to look for and a directory to look in
# must have a unique id - will collect all results from any folders with id in their names
collate.csv <- function(directory,file.name = "testStatistics",file_id,run_id,output_path){
  # Collect all the folders within the directory
  folder_paths <- list.dirs(directory, full.names = FALSE, recursive = FALSE)
  # Now reduce that to only get folders for the particular id of interest
  id_inds <- grep(file_id,folder_paths)
  id_folder_paths <- folder_paths[id_inds]
  # Make the path to each csv file using the file name
  csv_paths <- paste0(directory, id_folder_paths, "/",file.name,".csv")
  csv_paths_exist <- csv_paths[file.exists(csv_paths)] # only take paths that have an output
  # Set the number of rows in the dataframe (will equal the number of folders - one folder per simulation)
  num_rows <- length(csv_paths_exist)
  # Open all the csv files, store the results as a list
  output_list <- lapply(csv_paths_exist, read.csv, stringsAsFactors = FALSE)
  # Reduce the dataframe from a list into a matrix
  output_df <- Reduce(rbind, output_list)
  # Remove the "X" column (row number for smaller csv files)
  output_df <- output_df[,2:ncol(output_df)]
  # save the output dataframe
  write.csv(output_df,file=paste0(output_path,file_id,"_",file.name,"_collatedSimulationData_",run_id,".csv"))
  
  # now save the paths that don't have an output
  csv_paths_missing <- csv_paths[file.exists(csv_paths)==FALSE]
  write.csv(csv_paths_missing,file=paste0(output_path,file_id,"_",file.name,"_missingSimulations_filePaths_",run_id,".csv"))
}

collate.bootstraps <- function(directory, file.name, id, output.file.name){
  # Collect all the folders within the directory
  all_files <- list.files(directory, full.names = FALSE, recursive = TRUE)
  # Now reduce that to only get folders for the particular id of interest
  id_files <- all_files[grep(id,all_files)]
  # Get the files of interest and their full file paths
  csv_paths <- id_files[grep(file.name,id_files)]
  # Remove collated test statistic file from the list of csvs
  collated_name <- basename(output.file.name)
  csv_paths <- csv_paths[!grepl(collated_name,csv_paths)]
  if (length(csv_paths) > 0){
    csv_paths <- paste0(dirname(directory),"/",basename(directory),"/",csv_paths)
    # Set the number of rows in the dataframe (will equal the number of csv files: one per simulation + one for the original empirical alignment)
    num_rows <- length(csv_paths)
    # Open all the csv files, store the results as a list
    output_list <- lapply(csv_paths, read.csv, stringsAsFactors = FALSE)
    # Reduce the dataframe from a list into a matrix
    output_df <- Reduce(rbind, output_list)
    # save the output dataframe
    write.csv(output_df, file = output.file.name, row.names = FALSE)
    return(output_df)
  } else {
    return(NA)
  }
}

# Collect all the csv paths that have unfinished simulations
collate.missing.params <- function(directory,id){
  # Collect all the folders within the directory
  folder_paths <- list.files(directory)
  # Get the positions at which the simulations containing the id are at
  inds <- grep(id, folder_paths)
  # Get only the relevant folders (folders containing the id)
  csv_paths <- folder_paths[inds]
  # Make the path to each output csv file
  dir_paths <- paste0(directory,csv_paths)
  check_csv_paths <- paste0(directory, csv_paths, "/testStatistics.csv")
  params_paths <- dir_paths[file.exists(check_csv_paths)==FALSE] # only take paths that don't have a test statistic output
  # Set the number of rows in the dataframe (will equal the number of folders - one folder per simulation)
  num_rows <- length(params_paths)
  # add the params.csv file so that you can open and store them all
  params_paths <- paste0(params_paths,"/params.csv")
  # Open all the csv files, store the results as a list
  output_list <- lapply(params_paths, read.csv, stringsAsFactors = FALSE) # open all the parameter csv paths
  # Reduce the dataframe from a list into a matrix
  output_df <- Reduce(rbind, output_list)
  # Remove the "X" column (row number for smaller csv files)
  output_df <- output_df[,4:33]
  return(output_df)
}

# Given a folder containing a parameters file, extracts the run_id that was used when making that parameters file
extract.run.id <- function(folder){
  fs <- list.files(folder) # get names of all files in the folder
  f <- fs[grep("parameters",fs)][1] # take the first parameters file from the results folder
  f_split <- strsplit(f, split = "_")[[1]] # split the filename and take the last piece (the run id and the filetype)
  run_id <- gsub(".csv","",f_split[length(f_split)]) # replace filetype with "" to get run id
  return(run_id)
}


