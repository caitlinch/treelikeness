# R code to import and collate test statistic results, and to process the results

# Specify which file paths to use
run_location = "mac"
# run_location = "soma"

if (run_location == "mac"){
  # Set file paths etc
  raw_data_folder <- "/Users/caitlincherryh/Documents/Results/Output/"
  output_folder <- "/Users/caitlincherryh/Documents/Results/collatedOutput/"
  
  # Set working directory
  maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # for work computer
} else if (run_location == "soma") {
  # Set file paths etc
  raw_data_folder <- "/data/caitlin/treelikeness/output/"
  output_folder <- "/data/caitlin/treelikeness/results/"
  
  # Set working directory
  maindir <- "/data/caitlin/treelikeness/"
}

# Source files for functions
source(paste0(maindir,"code/func_process_data.R"))

# load required libraries
library(reshape2)

# Collate data for the four plots/sets of simulations and output each collated dataframe as a csv file
collate.csv(directory = raw_data_folder, file.name = "testStatistics", id = "plot1", output_path = output_folder)
collate.csv(directory = raw_data_folder, file.name = "testStatistics", id = "plot2", output_path = output_folder)
collate.csv(directory = raw_data_folder, file.name = "testStatistics", id = "plot3", output_path = output_folder)
plot4_ids <- paste0("plot4tree",1:9)
for (i in plot4_ids){
  collate.csv(directory = raw_data_folder, file.name = "testStatistics", id = i, output_path = output_folder)
  collate.csv(directory = raw_data_folder, file.name = "p_value", id = i, output_path = output_folder)
}

# Collate the files for plot4
csvs <- list.files(output_folder)
csvs <- csvs[grep("plot4tree",csvs)]
ts <- csvs[grep("testStatistics_collatedSimulationData",csvs)]
ts <- paste0(output_folder,ts)
ts_op_list <- lapply(ts, read.csv, stringsAsFactors = FALSE)
ts_op_df <- Reduce(rbind, ts_op_list)
ts_op_df <- ts_op_df[,2:ncol(ts_op_df)]
write.csv(ts_op_df,file=paste0(output_folder,"plot4_testStatistics_collatedSimulationData.csv"))
ps <- csvs[grep("p_value_collatedSimulationData",csvs)]
ps <- paste0(output_folder,ps)
ps_op_list <- lapply(ps, read.csv, stringsAsFactors = FALSE)
ps_op_df <- Reduce(rbind, ps_op_list)
ps_op_df <- ps_op_df[,2:ncol(ps_op_df)]
write.csv(ps_op_df,file=paste0(output_folder,"plot4_p_value_collatedSimulationData.csv"))

# Calculate the proportion of recombinant triplets and add it onto each set of simulations
id <- c("plot1_","plot2_","plot3_","plot4_")
csvs <- list.files(output_folder)
inds <- lapply(id,grep,csvs)
csvs <- csvs[unlist(inds)]
csvs <- csvs[grep("testStatistics_collatedSimulationData",csvs)]
csvs <- paste0(output_folder,csvs)
for (csv in csvs){
  df <- read.csv(csv, stringsAsFactors = FALSE)
  # divide the number of recombinant triplets detected by 3seq by the number of triplets tested
  # To find # of triplets: "In a set of 10 sequences, there are 720 unique parent–parent–child arrangements" - Boni et al (2007)
  # In other words: 6*choose(10,3) == 720
  # number of triplets tested will be 6* n choose k (if have a,b,c: a and b can be parents, b and c can be parents and a and c can be parents BUT each parent can be either P or Q)
  df["num_3seq_triplets"] <- 6 * choose(df$n_taxa, 3)
  df["proportion_recombinant_triplets"] <- df$X3SEQ_num_recombinant_triplets / df$num_3seq_triplets
  write.csv(df, file = csv)
}

# Reshape the data into long format
for (csv in csvs){
  df <- read.csv(csv, stringsAsFactors = FALSE)
  id_vars <- c("method","n_taxa","n_sites","tree_age","tree1","proportion_tree1","tree2","proportion_tree2","id")
  measure_vars <- c("PHI_observed","prop_resolved_quartets","proportion_recombinant_triplets","splittable_percentage","pdm_difference","neighbour_net")
  melt_df <- melt(df, id = id_vars, measure.vars = measure_vars)
  output_name <- gsub(".csv","_melted.csv",csv)
  write.csv(melt_df, file = output_name)
}
