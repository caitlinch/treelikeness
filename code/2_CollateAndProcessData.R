# R code to import and collate test statistic results, and to process the results

# Set file paths etc
raw_data_folder <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments2/"
output_folder <- "/Users/caitlincherryh/Documents/TestAlignmentResults/0_prelim_mk2/"

# Set working directory
maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # for work computer
setwd(maindir)

# Source files for functions
source(paste0(maindir,"code/func_process_data.R"))

# load required libraries
library(ggplot2)
library(reshape2)
library(gridExtra)

# Collate data for the three preliminary sets of simulations and output each collated dataframe as a csv file
external_df <- collate.csv(raw_data_folder, "external")
external_df <- simplify.SimBac(external_df)
file_name <- paste0(output_folder, "0_prelimMk2_collated_external.csv")
write.csv(external_df, file = file_name, row.names = FALSE)

internal_df <- collate.csv(raw_data_folder, "internal")
internal_df <- simplify.SimBac(internal_df)
file_name <- paste0(output_folder, "0_prelimMk2_collated_internal.csv")
write.csv(internal_df, file = file_name, row.names = FALSE)

phylo_df <- collate.csv(raw_data_folder, "pattern")
phylo_df <- simplify.phylo(phylo_df)
file_name <- paste0(output_folder, "0_prelimMk2_collated_pattern.csv")
write.csv(phylo_df, file = file_name, row.names = FALSE)

# divide the number of recombinant triplets detected by 3seq by the number of triplets tested
# number of triplets tested will be 6* n choose k (if have a,b,c: a and b can be parents, b and c can be parents and a and c can be parents BUT each parent can be either P or Q)
external_df["num_3seq_triplets"] <- 6 * choose(external_df$n_taxa, 3)
external_df["proportion_recombinant_triplets"] <- external_df$X3SEQ_num_recombinant_triplets / external_df$num_3seq_triplets

internal_df["num_3seq_triplets"] <- 6 * choose(internal_df$n_taxa, 3)
internal_df["proportion_recombinant_triplets"] <- internal_df$X3SEQ_num_recombinant_triplets / internal_df$num_3seq_triplets

phylo_df["num_3seq_triplets"] <- 6 * choose(phylo_df$n_taxa, 3)
phylo_df["proportion_recombinant_triplets"] <- phylo_df$X3SEQ_num_recombinant_triplets / phylo_df$num_3seq_triplets
phylo_df["proportion_tree2"] <- 1 - phylo_df$proportion_tree1 # get the proportion of the second tree (increasing proportion of recombination)

# Get only the relevant columns to plot
# PHI (pairwise homoplasy [convergence] index) - gives p-value of observing the sequences under the null hypothesis of no recombination calculated using the PHI statistic.
# Extract the observed PHI column and the other columns you want (including the new "proportion_recombinant_triplets" column)
external_pruned_df <- external_df[ , c("external_recombination","PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","pdm_average","split_decomposition","neighbour_net")]
internal_pruned_df <- internal_df[ , c("internal_recombination","PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","pdm_average","split_decomposition","neighbour_net")]
phylo_pruned_df <- phylo_df[ , c("proportion_tree2","PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","pdm_average","split_decomposition","neighbour_net")]

# Reshape the data into long format
plot_ex_df <- melt(external_pruned_df, id = c("external_recombination"))
plot_in_df <- melt(internal_pruned_df, id = c("internal_recombination"))
plot_phy_df <- melt(phylo_pruned_df, id = c("proportion_tree2"))

# Plot the data
# Plot of all stats on the same set of axes
ex_plot <- ggplot(plot_ex_df, aes(x = external_recombination, y = value, col = variable, shape = variable)) + geom_point(size = 2) + 
            labs(x = "External recombination (%)", y = "Statistic value", title = "All Test Statistics") + 
            stat_summary(fun.y = mean, geom="line", lwd=1, aes(col = variable)) +
            scale_colour_manual(name  = "Statistic", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                                labels = c("PHI", "Proportion Recombinant Triplets", "Proportion Resolved Quartets", 
                                           "TS1 - Splittable Percentage", "TS2a - Difference", "TS2b - Normalised Difference Mean", 
                                           "TS3a - Split Decomposition", "TS3b - NeighbourNet")) + 
            scale_shape_manual(name  = "Statistic", values = c(16, 17, 18, 15, 12, 8, 4, 1), 
                               labels = c("PHI", "Proportion Recombinant Triplets", "Proportion Resolved Quartets", 
                                          "TS1 - Splittable Percentage", "TS2a - Difference", "TS2b - Normalised Difference Mean", 
                                          "TS3a - Split Decomposition", "TS3b - NeighbourNet")) + 
            scale_linetype_manual(name  = "Statistic", values = c(1, 2, 3, 4, 5, 6, 1, 2), 
                                  labels = c("PHI", "Proportion Recombinant Triplets", "Proportion Resolved Quartets", 
                                  "TS1 - Splittable Percentage", "TS2a - Difference", "TS2b - Normalised Difference Mean", 
                                  "TS3a - Split Decomposition", "TS3b - NeighbourNet")) + theme_light()
plot_title <- paste0(output_folder, "externalRecombination_allTests.pdf")
ggsave(filename = plot_title, plot = ex_plot, dev = "pdf")
plot_title <- paste0(output_folder, "externalRecombination_allTests.png")
ggsave(filename = plot_title, plot = ex_plot, dev = "png")

# Plot each statistic on its own plot
ex_plot1<- ggplot(external_pruned_df, aes(x = external_recombination, y = PHI_observed)) + geom_point(size = 2, color = "#999999", shape = 16) + 
      labs(x = "External recombination (%)", y = "Statistic value", title = "Pairwise Homoplasy Index (PHI)") +
      stat_summary(fun.y = mean, geom="line", lwd=1, col = "#999999", lty = 1) + theme_light()

ex_plot2 <- ggplot(external_pruned_df, aes(x = external_recombination, y =proportion_recombinant_triplets)) + geom_point(size = 2, color = "#E69F00", shape = 17) + 
  labs(x = "External recombination (%)", y = "Statistic value", title = "Proportion Recombinant Triplets") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#E69F00", lty = 1) + theme_light()

ex_plot3 <- ggplot(external_pruned_df, aes(x = external_recombination, y =prop_resolved_quartets)) + geom_point(size = 2, color = "#56B4E9", shape = 18) + 
  labs(x = "External recombination (%)", y = "Statistic value", title = "Proportion of Resolved Quartets") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#56B4E9", lty = 1) + theme_light()

ex_plot4 <- ggplot(external_pruned_df, aes(x = external_recombination, y = splittable_percentage)) + geom_point(size = 2, color = "#009E73", shape = 15) + 
  labs(x = "External recombination (%)", y = "Statistic value", title = "TS1 - Splittable Percentage") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#009E73", lty = 1) + theme_light()

ex_plot5 <- ggplot(external_pruned_df, aes(x = external_recombination, y = pdm_difference)) + geom_point(size = 2, color = "#F0E442", shape = 12) + 
  labs(x = "External recombination (%)", y = "Statistic value", title = "TS2a - Difference") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#F0E442", lty = 1) + theme_light()

ex_plot6 <- ggplot(external_pruned_df, aes(x = external_recombination, y = pdm_average)) + geom_point(size = 2, color = "#0072B2", shape = 8) + 
  labs(x = "External recombination (%)", y = "Statistic value", title = "TS2b - Normalised Difference Mean") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#0072B2", lty = 1) + theme_light()

ex_plot7 <- ggplot(external_pruned_df, aes(x = external_recombination, y = split_decomposition)) + geom_point(size = 2, color = "#D55E00", shape = 4) + 
  labs(x = "External recombination (%)", y = "Statistic value", title = "TS3a - Split Decomposition") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#D55E00", lty = 1) + theme_light()

ex_plot8 <- ggplot(external_pruned_df, aes(x = external_recombination, y = neighbour_net)) + geom_point(size = 2, color = "#CC79A7", shape = 1) + 
  labs(x = "External recombination (%)", y = "Statistic value", title = "TS3b - NeighbourNet") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#CC79A7", lty = 1) + theme_light()

# Assemble all the plots together and save them
combined_plot <- grid.arrange(ex_plot1, ex_plot2, ex_plot3, ex_plot4, ex_plot5, ex_plot6, ex_plot7, ex_plot8 , nrow = 2, ncol = 4)
plot_title <- paste0(output_folder, "externalRecombination_MosaicPlots.pdf")
ggsave(filename = plot_title, plot = combined_plot, dev = "pdf", width = 32.44, height = 14.88, units = "in")
plot_title <- paste0(output_folder, "externalRecombination_MosaicPlots.png")
ggsave(filename = plot_title, plot = combined_plot, dev = "png", width = 32.44, height = 14.88, units = "in")


# Repeat this for the internal recombination plots
in_plot <- ggplot(plot_in_df, aes(x = internal_recombination, y = value, col = variable, shape = variable)) + geom_point(size = 2) + 
  labs(x = "Internal recombination (%)", y = "Statistic value", title = "All Test Statistics") + 
  stat_summary(fun.y = mean, geom="line", lwd=1, aes(col = variable)) +
  scale_colour_manual(name  = "Statistic", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                      labels = c("PHI", "Proportion Recombinant Triplets", "Proportion Resolved Quartets", 
                                 "TS1 - Splittable Percentage", "TS2a - Difference", "TS2b - Normalised Difference Mean", 
                                 "TS3a - Split Decomposition", "TS3b - NeighbourNet")) + 
  scale_shape_manual(name  = "Statistic", values = c(16, 17, 18, 15, 12, 8, 4, 1), 
                     labels = c("PHI", "Proportion Recombinant Triplets", "Proportion Resolved Quartets", 
                                "TS1 - Splittable Percentage", "TS2a - Difference", "TS2b - Normalised Difference Mean", 
                                "TS3a - Split Decomposition", "TS3b - NeighbourNet")) + 
  scale_linetype_manual(name  = "Statistic", values = c(1, 2, 3, 4, 5, 6, 1, 2), 
                        labels = c("PHI", "Proportion Recombinant Triplets", "Proportion Resolved Quartets", 
                                   "TS1 - Splittable Percentage", "TS2a - Difference", "TS2b - Normalised Difference Mean", 
                                   "TS3a - Split Decomposition", "TS3b - NeighbourNet")) + theme_light()
plot_title <- paste0(output_folder, "internalRecombination_allTests.pdf")
ggsave(filename = plot_title, plot = in_plot, dev = "pdf")
plot_title <- paste0(output_folder, "internalRecombination_allTests.png")
ggsave(filename = plot_title, plot = in_plot, dev = "png")

# Plot each statistic on its own plot

in_plot1<- ggplot(internal_pruned_df, aes(x = internal_recombination, y = PHI_observed)) + geom_point(size = 2, color = "#999999", shape = 16) + 
  labs(x = "Internal recombination (%)", y = "Statistic value", title = "Pairwise Homoplasy Index (PHI)") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#999999", lty = 1) + theme_light()

in_plot2 <- ggplot(internal_pruned_df, aes(x = internal_recombination, y = proportion_recombinant_triplets)) + geom_point(size = 2, color = "#E69F00", shape = 17) + 
  labs(x = "Internal recombination (%)", y = "Statistic value", title = "Proportion Recombinant Triplets") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#E69F00", lty = 1) + theme_light()

in_plot3 <- ggplot(internal_pruned_df, aes(x = internal_recombination, y = prop_resolved_quartets)) + geom_point(size = 2, color = "#56B4E9", shape = 18) + 
  labs(x = "Internal recombination (%)", y = "Statistic value", title = "Proportion of Resolved Quartets") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#56B4E9", lty = 1) + theme_light()

in_plot4 <- ggplot(internal_pruned_df, aes(x = internal_recombination, y = splittable_percentage)) + geom_point(size = 2, color = "#009E73", shape = 15) + 
  labs(x = "Internal recombination (%)", y = "Statistic value", title = "TS1 - Splittable Percentage") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#009E73", lty = 1) + theme_light()

in_plot5 <- ggplot(internal_pruned_df, aes(x = internal_recombination, y = pdm_difference)) + geom_point(size = 2, color = "#F0E442", shape = 12) + 
  labs(x = "Internal recombination (%)", y = "Statistic value", title = "TS2a - Difference") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#F0E442", lty = 1) + theme_light()

in_plot6 <- ggplot(internal_pruned_df, aes(x = internal_recombination, y = pdm_average)) + geom_point(size = 2, color = "#0072B2", shape = 8) + 
  labs(x = "Internal recombination (%)", y = "Statistic value", title = "TS2b - Normalised Difference Mean") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#0072B2", lty = 1) + theme_light()

in_plot7 <- ggplot(internal_pruned_df, aes(x = internal_recombination, y = split_decomposition)) + geom_point(size = 2, color = "#D55E00", shape = 4) + 
  labs(x = "Internal recombination (%)", y = "Statistic value", title = "TS3a - Split Decomposition") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#D55E00", lty = 1) + theme_light()

in_plot8 <- ggplot(internal_pruned_df, aes(x = internal_recombination, y = neighbour_net)) + geom_point(size = 2, color = "#CC79A7", shape = 1) + 
  labs(x = "Internal recombination (%)", y = "Statistic value", title = "TS3b - NeighbourNet") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#CC79A7", lty = 1) + theme_light()

# Assemble all the plots together and save them
combined_plot <- grid.arrange(in_plot1, in_plot2, in_plot3, in_plot4, in_plot5, in_plot6, in_plot7, in_plot8 , nrow = 2, ncol = 4)
plot_title <- paste0(output_folder, "internalRecombination_MosaicPlots.pdf")
ggsave(filename = plot_title, plot = combined_plot, dev = "pdf", width = 32.44, height = 14.88, units = "in")
plot_title <- paste0(output_folder, "internalRecombination_MosaicPlots.png")
ggsave(filename = plot_title, plot = combined_plot, dev = "png", width = 32.44, height = 14.88, units = "in")

# Repeat for the phylo plots
phy_plot <- ggplot(plot_phy_df, aes(x = proportion_tree2, y = value, col = variable, shape = variable)) + geom_point(size = 2) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "All Test Statistics") + 
  stat_summary(fun.y = mean, geom="line", lwd=1, aes(col = variable)) +
  scale_colour_manual(name  = "Statistic", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                      labels = c("PHI", "Proportion Recombinant Triplets", "Proportion Resolved Quartets", 
                                 "TS1 - Splittable Percentage", "TS2a - Difference", "TS2b - Normalised Difference Mean", 
                                 "TS3a - Split Decomposition", "TS3b - NeighbourNet")) + 
  scale_shape_manual(name  = "Statistic", values = c(16, 17, 18, 15, 12, 8, 4, 1), 
                     labels = c("PHI", "Proportion Recombinant Triplets", "Proportion Resolved Quartets", 
                                "TS1 - Splittable Percentage", "TS2a - Difference", "TS2b - Normalised Difference Mean", 
                                "TS3a - Split Decomposition", "TS3b - NeighbourNet")) + 
  scale_linetype_manual(name  = "Statistic", values = c(1, 2, 3, 4, 5, 6, 1, 2), 
                        labels = c("PHI", "Proportion Recombinant Triplets", "Proportion Resolved Quartets", 
                                   "TS1 - Splittable Percentage", "TS2a - Difference", "TS2b - Normalised Difference Mean", 
                                   "TS3a - Split Decomposition", "TS3b - NeighbourNet")) + theme_light()
plot_title <- paste0(output_folder, "phylo_allTests_mean.pdf")
ggsave(filename = plot_title, plot = phy_plot, dev = "pdf")
plot_title <- paste0(output_folder, "phylo_allTests_mean.png")
ggsave(filename = plot_title, plot = phy_plot, dev = "png")

# Plot each statistic on its own plot

phy_plot1<- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y = PHI_observed)) + geom_point(size = 2, color = "#999999", shape = 16) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "Pairwise Homoplasy Index (PHI)") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#999999", lty = 1) + theme_light()

phy_plot2 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y =proportion_recombinant_triplets)) + geom_point(size = 2, color = "#E69F00", shape = 17) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "Proportion Recombinant Triplets") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#E69F00", lty = 1) + theme_light()

phy_plot3 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y =prop_resolved_quartets)) + geom_point(size = 2, color = "#56B4E9", shape = 18) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "Proportion of Resolved Quartets") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#56B4E9", lty = 1) + theme_light()

phy_plot4 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y = splittable_percentage)) + geom_point(size = 2, color = "#009E73", shape = 15) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "TS1 - Splittable Percentage") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#009E73", lty = 1) + theme_light()

phy_plot5 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y = pdm_difference)) + geom_point(size = 2, color = "#F0E442", shape = 12) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "TS2a - Difference") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#F0E442", lty = 1) + theme_light()

phy_plot6 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y = pdm_average)) + geom_point(size = 2, color = "#0072B2", shape = 8) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "TS2b - Normalised Difference Mean") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#0072B2", lty = 1) + theme_light()

phy_plot7 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y = split_decomposition)) + geom_point(size = 2, color = "#D55E00", shape = 4) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "TS3a - Split Decomposition") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#D55E00", lty = 1) + theme_light()

phy_plot8 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y = neighbour_net)) + geom_point(size = 2, color = "#CC79A7", shape = 1) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "TS3b - NeighbourNet") +
  stat_summary(fun.y = mean, geom="line", lwd=1, col = "#CC79A7", lty = 1) + theme_light()

# Assemble all the plots together and save them
combined_plot <- grid.arrange(phy_plot1, phy_plot2, phy_plot3, phy_plot4, phy_plot5, phy_plot6, phy_plot7, phy_plot8 , nrow = 2, ncol = 4)
plot_title <- paste0(output_folder, "phylo_MosaicPlots_mean.pdf")
ggsave(filename = plot_title, plot = combined_plot, dev = "pdf", width = 32.44, height = 14.88, units = "in")
plot_title <- paste0(output_folder, "phylo_MosaicPlots_mean.png")
ggsave(filename = plot_title, plot = combined_plot, dev = "png", width = 32.44, height = 14.88, units = "in")

# Repeat for the phylo plots but with a linear model
phy_plot <- ggplot(plot_phy_df, aes(x = proportion_tree2, y = value, col = variable, shape = variable)) + geom_point(size = 2) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "All Test Statistics") + 
  geom_smooth(method = "lm", aes(col = variable)) +
  scale_colour_manual(name  = "Statistic", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                      labels = c("PHI", "Proportion Recombinant Triplets", "Proportion Resolved Quartets", 
                                 "TS1 - Splittable Percentage", "TS2a - Difference", "TS2b - Normalised Difference Mean", 
                                 "TS3a - Split Decomposition", "TS3b - NeighbourNet")) + 
  scale_shape_manual(name  = "Statistic", values = c(16, 17, 18, 15, 12, 8, 4, 1), 
                     labels = c("PHI", "Proportion Recombinant Triplets", "Proportion Resolved Quartets", 
                                "TS1 - Splittable Percentage", "TS2a - Difference", "TS2b - Normalised Difference Mean", 
                                "TS3a - Split Decomposition", "TS3b - NeighbourNet")) + 
  scale_linetype_manual(name  = "Statistic", values = c(1, 2, 3, 4, 5, 6, 1, 2), 
                        labels = c("PHI", "Proportion Recombinant Triplets", "Proportion Resolved Quartets", 
                                   "TS1 - Splittable Percentage", "TS2a - Difference", "TS2b - Normalised Difference Mean", 
                                   "TS3a - Split Decomposition", "TS3b - NeighbourNet")) + theme_light()
plot_title <- paste0(output_folder, "phylo_allTests_lm.pdf")
ggsave(filename = plot_title, plot = phy_plot, dev = "pdf")
plot_title <- paste0(output_folder, "phylo_allTests_lm.png")
ggsave(filename = plot_title, plot = phy_plot, dev = "png")

# Plot each statistic on its own plot
phy_plot1<- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y = PHI_observed)) + geom_point(size = 2, color = "#999999", shape = 16) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "Pairwise Homoplasy Index (PHI)") +
  geom_smooth(method = "lm", aes(col = "#999999")) + theme_light()

phy_plot2 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y =proportion_recombinant_triplets)) + geom_point(size = 2, color = "#E69F00", shape = 17) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "Proportion Recombinant Triplets") +
  geom_smooth(method = "lm", aes(col = "#E69F00")) + theme_light()

phy_plot3 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y =prop_resolved_quartets)) + geom_point(size = 2, color = "#56B4E9", shape = 18) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "Proportion of Resolved Quartets") +
  geom_smooth(method = "lm", aes(col = "#56B4E9")) + theme_light()

phy_plot4 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y = splittable_percentage)) + geom_point(size = 2, color = "#009E73", shape = 15) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "TS1 - Splittable Percentage") +
  geom_smooth(method = "lm", aes(col = "#009E73")) + theme_light()

phy_plot5 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y = pdm_difference)) + geom_point(size = 2, color = "#F0E442", shape = 12) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "TS2a - Difference") +
  geom_smooth(method = "lm", aes(col = "#F0E442")) + theme_light()

phy_plot6 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y = pdm_average)) + geom_point(size = 2, color = "#0072B2", shape = 8) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "TS2b - Normalised Difference Mean") +
  geom_smooth(method = "lm", aes(col = "#0072B2")) + theme_light()

phy_plot7 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y = split_decomposition)) + geom_point(size = 2, color = "#D55E00", shape = 4) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "TS3a - Split Decomposition") +
  geom_smooth(method = "lm", aes(col = "#D55E00")) + theme_light()

phy_plot8 <- ggplot(phylo_pruned_df, aes(x = proportion_tree2, y = neighbour_net)) + geom_point(size = 2, color = "#CC79A7", shape = 1) + 
  labs(x = "Proportion Tree 2 (%)", y = "Statistic value", title = "TS3b - NeighbourNet") +
  geom_smooth(method = "lm", aes(col = "#CC79A7")) + theme_light()

# Assemble all the plots together and save them
combined_plot <- grid.arrange(phy_plot1, phy_plot2, phy_plot3, phy_plot4, phy_plot5, phy_plot6, phy_plot7, phy_plot8 , nrow = 2, ncol = 4)
plot_title <- paste0(output_folder, "phylo_MosaicPlots_lm.pdf")
ggsave(filename = plot_title, plot = combined_plot, dev = "pdf", width = 32.44, height = 14.88, units = "in")
plot_title <- paste0(output_folder, "phylo_MosaicPlots_lm.png")
ggsave(filename = plot_title, plot = combined_plot, dev = "png", width = 32.44, height = 14.88, units = "in")

