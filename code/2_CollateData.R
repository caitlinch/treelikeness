# R code to import and collate test statistic results, and to process the results
# Sourcing this file will extract, collect, format the results from the simulations, and perform some additional calculations
# Final result is four dataframes (three containing experimental results and one containing results from the parametric bootstrap) 
# that can be used for analysis



##### Step 1: Open packages #####
library(reshape2)



##### Step 2: Uncomment and set the file paths for output folders, executables, and the identifying name for this run #####
# op_folder <- the folder where simulated alignments and output from analysis (e.g. IQ-Tree output files, 3seq output files, test statistic csvs) 
#              are placed. MUST be the same folder as in Part 1, as it looks for these files to extract test statistics and other information.
# results_folder <- the folder where the result csvs will be placed (I use same results_folder in Parts 1-4)
# maindir <- "treelikeness" repository location
# run_id <- if "run.id  = FALSE", program extracts run_id from input parameter file names 
#        <- otherwise, run_id will be set to whatever the user inputs here (e.g. "run_id = 'replicateAnalysis' ")

#__________________________________________File paths/parameters: update these for your own machine______________________________________
# op_folder <- ""
# results_folder <- ""
# maindir <- ""
# run_id <- extract.run.id(results_folder) 
#____________________________________________________________________________________________________________________________________

#__________________________________________Caitlin's paths (delete these if you're not Caitlin)______________________________________
op_folder <- "/data/caitlin/treelikeness/output_20200304/"
results_folder <- "/data/caitlin/treelikeness/results_20200304/"
maindir <- "/data/caitlin/treelikeness/"
run_id = FALSE
#____________________________________________________________________________________________________________________________________




##### Step 3: Source function files #####
source(paste0(maindir,"code/func_process_data.R"))
# Extract run.id from the results folder name so the whole analysis has the same run.id
if (run_id == "FALSE"){
  run_id <- extract.run.id(results_folder)
}



##### Step 4: Collect test statistics from output, collate into a single file and write dataframes #####
# Collate data for the three sets of simulations and output each collated dataframe as a csv file
collate.csv(directory = op_folder, file.name = "testStatistics", file_id = "exp2", run_id = run_id, output_path = results_folder)
collate.csv(directory = op_folder, file.name = "p_value", file_id = "exp2", run_id = run_id, output_path = results_folder)
collate.csv(directory = op_folder, file.name = "testStatistics", file_id = "exp3", run_id = run_id, output_path = results_folder)
collate.csv(directory = op_folder, file.name = "p_value", file_id = "exp3", run_id = run_id, output_path = results_folder)



##### Step 5: Calculate additional test statistics and format dataframes #####
# Will result in one df containing test statistic values for each of 3 experiments, and one df containing p-value results for the third experiment
# Extract the collated file names 
id <- c("exp2_","exp3_")
csvs <- list.files(results_folder)
inds <- lapply(id,grep,csvs)
csvs <- csvs[unlist(inds)]
csvs <- c(csvs[grep("testStatistics_collatedSimulationData",csvs)], csvs[grep("p_value_collatedSimulationData",csvs)])
csvs <- csvs[grep("melted",csvs,invert=TRUE)] # If this is a rerun, take only the collated data files and ignore the "melted" (final result of Part 2) files
csvs <- paste0(results_folder,csvs)
# Add columns for calculating 3seq-based test statistic (for exp dfs only) and to detail information about the recombination event (all dfs)
for (csv in csvs){
  df <- read.csv(csv, stringsAsFactors = FALSE)
  # Calculate the number of triplets tested by 3seq and use to calculate the proportion of recombinant triplets for the experiment dfs
  if ("testStatistics" %in% strsplit(csv,"_")[[1]]){
    # divide the number of recombinant triplets detected by 3seq by the number of triplets tested
    # To find # of triplets: "In a set of 10 sequences, there are 720 unique parent–parent–child arrangements" - Boni et al (2007)
    # In other words: 6*choose(10,3) == 720
    # number of triplets tested will be 6* n choose k (if have a,b,c: a and b can be parents, b and c can be parents and a and c can be parents BUT each parent can be either P or Q)
    df["X3SEQ_num_3seq_triplets"] <- 6 * choose(df$n_taxa, 3)
    df["X3SEQ_proportion_recombinant_triplets"] <- df$X3SEQ_num_recombinant_triplets / df$num_3seq_triplets
  }
  # Extract information about the type of recombination event from the tree names and add to columns so it can be used for analysis (for all dfs)
  vector <- 1:length(df$tree2)
  close_inds <- grep("close",df$tree2)
  divergent_inds <- grep("divergent",df$tree2)
  ancient_inds <- grep("ancient",df$tree2)
  noevent_inds <- grep("LHS",df$tree2)
  tree2_event_position <- vector
  tree2_event_position[close_inds] <- "close"
  tree2_event_position[divergent_inds] <- "divergent"
  tree2_event_position[ancient_inds] <- "ancient"
  tree2_event_position[noevent_inds] <- "none"
  balanced_inds <- grep("balanced",df$tree1)
  intermediate_inds <- grep("intermediate",df$tree1)
  unbalanced_inds <- grep("unbalanced",df$tree1)
  tree1_shape <- vector
  tree1_shape[balanced_inds] <- "balanced"
  tree1_shape[intermediate_inds] <- "intermediate"
  tree1_shape[unbalanced_inds] <- "unbalanced"
  balanced_inds <- grep("balanced",df$tree2)
  intermediate_inds <- grep("intermediate",df$tree2)
  unbalanced_inds <- grep("unbalanced",df$tree2)
  tree2_shape <- vector
  tree2_shape[balanced_inds] <- "balanced"
  tree2_shape[intermediate_inds] <- "intermediate"
  tree2_shape[unbalanced_inds] <- "unbalanced"
  nonreciprocal_inds <- grep("_nonreciprocal",df$tree2)
  noevent_inds <- grep("LHS",df$tree2)
  reciprocal_inds <- grep("_reciprocal",df$tree2)
  tree2_event_type <- vector
  tree2_event_type[nonreciprocal_inds] <- "nonreciprocal"
  tree2_event_type[noevent_inds] <- "none"
  tree2_event_type[reciprocal_inds] <- "reciprocal"
  tree2_eventnum <- vector
  for (i in 1:8){
    temp_inds <- grep(paste0(i,"event"),df$tree2)
    tree2_eventnum[temp_inds] <- i
  }
  noevent_inds <- grep("LHS",df$tree2)
  tree2_eventnum[noevent_inds] <- 0
  tree2_eventnum <- as.numeric(tree2_eventnum)
  df["tree2_event_position"] <- tree2_event_position
  df["tree1_tree_shape"] <- tree1_shape
  df["tree2_tree_shape"] <- tree2_shape
  df["tree2_event_type"] <- tree2_event_type
  df["number_of_events"] <- tree2_eventnum
  write.csv(df, file = csv, row.names = FALSE)
}



##### Step 6: Reshape the data into long format and write dataframes #####
# Reshape experiment dfs into melted (long) format
ts_csvs <- csvs[grep("testStatistics",csvs)]
for (csv in ts_csvs){
  df <- read.csv(csv, stringsAsFactors = FALSE)
  id_vars <- c("n_taxa","n_sites","tree_age","tree1_tree_shape","proportion_tree1","tree2_event_position","tree2_event_type","tree2_tree_shape","proportion_tree2","number_of_events","id")
  measure_vars <- c("PHI_observed","X3SEQ_proportion_recombinant_triplets","LM_prop_resolved_quartets","tree_proportion","mean_delta_q", "median_delta_q")
  melt_df <- melt(df, id = id_vars, measure.vars = measure_vars)
  output_name <- gsub(".csv","_melted.csv",csv)
  write.csv(melt_df, file = output_name, row.names = FALSE)
}

# Reshape p value dfs into melted (long) format
p_csvs <- csvs[grep("p_value",csvs)]
for (csv in p_csvs){
  df <- read.csv(csv, stringsAsFactors = FALSE)
  id_vars <- c("n_taxa","n_sites","tree_age","tree1_tree_shape","proportion_tree1","tree2_event_position","tree2_event_type","tree2_tree_shape",
               "proportion_tree2","number_of_events","id")
  measure_vars <- c("PHI_p_value","X3Seq_p_value","LM_p_value", "tree_proportion_p_value","mean_delta_q_p_value", 
                    "median_delta_q_p_value")
  melt_df <- melt(df, id = id_vars, measure.vars = measure_vars)
  output_name <- gsub(".csv","_melted.csv",csv)
  write.csv(melt_df, file = output_name, row.names = FALSE)
}

