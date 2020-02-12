# R code to import and collate test statistic results, and to process the results
# Sourcing this file will extract, collect, format the results from the simulations, and perform some additional calculations

##### Step 1: Open packages #####
library(reshape2)



##### Step 2: Specify file paths #####
# op_folder <- the folder where simulated alignments and output from analysis (e.g. IQ-Tree output files, 3seq output files, test statistic csvs) 
#              are placed. MUST be the same folder as in Part 1, as it looks for these files to extract test statistics and other information.
# results_folder <- the folder where the result csvs will be placed (I use same results_folder in Parts 1 and 2)
# maindir <- "treelikeness" repository location

# op_folder <- ""
# results_folder <- ""
# maindir <- ""


run_location = "mac"
# run_location = "soma"

if (run_location == "mac"){
  # Set file paths etc
  op_folder <- "/Users/caitlincherryh/Documents/Honours/TestAlignmentResults/9_MStests/001_test/op/"
  results_folder <- "/Users/caitlincherryh/Documents/Honours/TestAlignmentResults/9_MStests/001_test/results/"
  maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # for work computer
} else if (run_location == "soma") {
  # Set file paths etc
  op_folder <- "/data/caitlin/treelikeness/output_20190411/"
  results_folder <- "/data/caitlin/treelikeness/results_20190411/"
  maindir <- "/data/caitlin/treelikeness/"
}



##### Step 3: Source function files #####
source(paste0(maindir,"code/func_process_data.R"))



##### Step 4: Collect test statistics from output and collate into a single file #####
# Collate data for the four plots/sets of simulations and output each collated dataframe as a csv file
collate.csv(directory = op_folder, file.name = "testStatistics", id = "exp1", output_path = results_folder)
collate.csv(directory = op_folder, file.name = "testStatistics", id = "exp2", output_path = results_folder)
collate.csv(directory = op_folder, file.name = "testStatistics", id = "exp3", output_path = results_folder)
collate.csv(directory = op_folder, file.name = "p_value", id = "exp3", output_path = results_folder)



##### Step 5: Calculate additional test statistics and format dataframes #####
# Calculate the proportion of recombinant triplets and add it onto each set of simulations
id <- c("exp1_","exp2_","exp3_")
csvs <- list.files(results_folder)
inds <- lapply(id,grep,csvs)
csvs <- csvs[unlist(inds)]
csvs <- c(csvs[grep("testStatistics_collatedSimulationData",csvs)], csvs[grep("p_value_collatedSimulationData",csvs)])
csvs <- paste0(results_folder,csvs)
for (csv in csvs[1:3]){
  df <- read.csv(csv, stringsAsFactors = FALSE)
  # divide the number of recombinant triplets detected by 3seq by the number of triplets tested
  # To find # of triplets: "In a set of 10 sequences, there are 720 unique parent–parent–child arrangements" - Boni et al (2007)
  # In other words: 6*choose(10,3) == 720
  # number of triplets tested will be 6* n choose k (if have a,b,c: a and b can be parents, b and c can be parents and a and c can be parents BUT each parent can be either P or Q)
  df["num_3seq_triplets"] <- 6 * choose(df$n_taxa, 3)
  df["proportion_recombinant_triplets"] <- df$X3SEQ_num_recombinant_triplets / df$num_3seq_triplets
  write.csv(df, file = csv, row.names = FALSE)
}

for (csv in csvs){
  df <- read.csv(csv, stringsAsFactors = FALSE)
  # divide the number of recombinant triplets detected by 3seq by the number of triplets tested
  # To find # of triplets: "In a set of 10 sequences, there are 720 unique parent–parent–child arrangements" - Boni et al (2007)
  # In other words: 6*choose(10,3) == 720
  # number of triplets tested will be 6* n choose k (if have a,b,c: a and b can be parents, b and c can be parents and a and c can be parents BUT each parent can be either P or Q)
  df["num_3seq_triplets"] <- 6 * choose(df$n_taxa, 3)
  df["proportion_recombinant_triplets"] <- df$X3SEQ_num_recombinant_triplets / df$num_3seq_triplets
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
for (csv in csvs[1:3]){
  df <- read.csv(csv, stringsAsFactors = FALSE)
  id_vars <- c("n_taxa","n_sites","tree_age","tree1_tree_shape","proportion_tree1","tree2_event_position","tree2_event_type","tree2_tree_shape","proportion_tree2","number_of_events","id")
  measure_vars <- c("PHI_observed","prop_resolved_quartets","proportion_recombinant_triplets","splittable_percentage","pdm_difference","neighbour_net_untrimmed","neighbour_net_trimmed",
                    "split_decomposition_untrimmed","split_decomposition_trimmed","mean_delta_q", "median_delta_q", "mode_delta_q")
  melt_df <- melt(df, id = id_vars, measure.vars = measure_vars)
  output_name <- gsub(".csv","_melted.csv",csv)
  write.csv(melt_df, file = output_name, row.names = FALSE)
}

# Reshape p value df into melted (long) format
df <- read.csv(csvs[4], stringsAsFactors = FALSE)
id_vars <- c("n_taxa","n_sites","tree_age","tree1_tree_shape","proportion_tree1","tree2_event_position","tree2_event_type","tree2_tree_shape","proportion_tree2","number_of_events","id")
measure_vars <- c("PHI_p_value","PHI_observed_p_value","X3Seq_p_value","num_recombinant_sequences_p_value","likelihood_mapping_p_value","splittable_percentage_p_value",
                  "pdm_difference_p_value","neighbour_net_untrimmed_p_value", "neighbour_net_trimmed_p_value","split_decomposition_untrimmed_p_value","split_decomposition_trimmed_p_value",
                  "mean_delta_q_p_value", "median_delta_q_p_value","mode_delta_q_p_value")
melt_df <- melt(df, id = id_vars, measure.vars = measure_vars)
output_name <- paste0(results_folder,"plot4_p_value_collatedSimulationData_melted.csv")
write.csv(melt_df, file = output_name, row.names = FALSE)
