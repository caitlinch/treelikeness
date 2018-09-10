# R code to import and collate test statistic results, and to process the results

# Specify which file paths to use
# run_location = "mac"
run_location = "soma"

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
library(ggplot2)
library(reshape2)
library(gridExtra)

# Collate data for the four plots/sets of simulations and output each collated dataframe as a csv file
collate.csv(directory = raw_data_folder, file.name = "testStatistics", id = "plot1", output_path = output_folder)
collate.csv(directory = raw_data_folder, file.name = "testStatistics", id = "plot2", output_path = output_folder)
collate.csv(directory = raw_data_folder, file.name = "testStatistics", id = "plot3", output_path = output_folder)
plot4_ids <- paste0("plot4tree",1:9)
for (i in plot4_ids){
  collate.csv(directory = raw_data_folder, file.name = "testStatistics", id = i, output_path = output_folder)
  collate.csv(directory = raw_data_folder, file.name = "p_value", id = i, output_path = output_folder)
}




# divide the number of recombinant triplets detected by 3seq by the number of triplets tested
# number of triplets tested will be 6* n choose k (if have a,b,c: a and b can be parents, b and c can be parents and a and c can be parents BUT each parent can be either P or Q)
# phylo_df["num_3seq_triplets"] <- 6 * choose(phylo_df$n_taxa, 3)
# phylo_df["proportion_recombinant_triplets"] <- phylo_df$X3SEQ_num_recombinant_triplets / phylo_df$num_3seq_triplets
# phylo_df["proportion_tree2"] <- 1 - phylo_df$proportion_tree1 # get the proportion of the second tree (increasing proportion of recombination)
 
# Get only the relevant columns to plot
# PHI (pairwise homoplasy [convergence] index) - gives p-value of observing the sequences under the null hypothesis of no recombination calculated using the PHI statistic.
# Extract the observed PHI column and the other columns you want (including the new "proportion_recombinant_triplets" column)
#external_pruned_df <- external_df[ , c("external_recombination","PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","pdm_average","split_decomposition","neighbour_net")]


# # Reshape the data into long format
# plot_ex_df <- melt(external_pruned_df, id = c("external_recombination"))
