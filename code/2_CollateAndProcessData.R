# R code to import and collate test statistic results, and to process the results

# Set file paths etc
output_folder <- "/Users/caitlincherryh/Documents/TestAlignmentResults/0_prelim_mk2/"

# Set working directory
maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # for work computer
setwd(maindir)

# Source files for functions
source(paste0(maindir,"code/func_process_data.R"))

# load required libraries
library(ggplot2)

# Collate data for the three preliminary sets of simulations and output each collated dataframe as a csv file
external_df <- collate.csv("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments2/","external")
external_df <- simplify.SimBac(external_df)
file_name <- paste0(output_folder,"0_prelimMk2_collated_external.csv")
write.csv(external_df, file = file_name, row.names = FALSE)

internal_df <- collate.csv("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments2/","internal")
internal_df <- simplify.SimBac(internal_df)
file_name <- paste0(output_folder,"0_prelimMk2_collated_internal.csv")
write.csv(internal_df, file = file_name, row.names = FALSE)

phylo_df <- collate.csv("/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments2/","pattern")
phylo_df <- simplify.phylo(phylo_df)
file_name <- paste0(output_folder,"0_prelimMk2_collated_pattern.csv")
write.csv(phylo_df, file = file_name, row.names = FALSE)

# divide the number of recombinant triplets detected by 3seq by the number of triplets tested
# number of triplets tested will be 6* n choose k (if have a,b,c: a and b can be parents, b and c can be parents and a and c can be parents BUT each parent can be either P or Q)
external_df["num_3seq_triplets"] <- 6*choose(external_df$n_taxa,3)
external_df["proportion_recombinant_triplets"] <- external_df$X3SEQ_num_recombinant_triplets/external_df$num_3seq_triplets

internal_df["num_3seq_triplets"] <- 6*choose(internal_df$n_taxa,3)
internal_df["proportion_recombinant_triplets"] <- internal_df$X3SEQ_num_recombinant_triplets/internal_df$num_3seq_triplets

phylo_df["num_3seq_triplets"] <- 6*choose(phylo_df$n_taxa,3)
phylo_df["proportion_recombinant_triplets"] <- phylo_df$X3SEQ_num_recombinant_triplets/phylo_df$num_3seq_triplets

# Get only the relevant columns to plot
# PHI (pairwise homoplasy [convergence] index) - gives p-value of observing the sequences under the null hypothesis of no recombination calculated using the PHI statistic.
# Extract the observed PHI column 
ex_plot_df <- external_df[,c("external_recombination","PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","pdm_average","split_decomposition","neighbour_net")]
in_plot_df <- internal_df[,c("internal_recombination","PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","pdm_average","split_decomposition","neighbour_net")]
phy_plot_df <- phylo_df[,c("proportion_tree1","PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","pdm_average","split_decomposition","neighbour_net")]






