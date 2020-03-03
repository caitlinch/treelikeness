# R code to import dataframes of test statistic and statistical test results and create plots to summarise the results
# Sourcing this file will open four dataframes and output a number of plots
# Final result is a number of plots displaying test statistic values under perturbation of various simulation factors



##### Step 1: Open packages #####
library(ggplot2)
library(gridExtra)
library(viridis)



##### Step 2: Specify file paths #####
# results_folder <- the folder where the result csvs will be placed (I use same results_folder in Parts 1 - 4.)
# plots_folder <- the folder where the plots will be stored
# maindir <- "treelikeness" repository location
# run_id <- program extracts run_id from input parameter file names

results_folder <- "/Users/caitlincherryh/Documents/Honours/Results/simulations_20190411/collatedOutput_2020/"
plots_folder <- "/Users/caitlincherryh/Documents/Honours/Results/simulations_20190411/plots_2020/"
maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/"
run_id <- extract.run.id(results_folder)



##### Step 3: Source function files #####
source(paste0(maindir,"code/func_process_data.R"))



##### Step 4: Open dataframes #####
# List all the files in the results_folder
melt_files <- list.files(results_folder)[grep("melted",list.files(results_folder))]
# Using the list of files in the results_folder, extract the melted csvs for each experiment (they are ready to plot!)
# p1:p3 are the dataframes for the test statistics calculated and estimated in experiments 1:3 respectively
p1_df <- read.csv(paste0(results_folder,melt_files[grep("exp1_testStatistics",melt_files)]), stringsAsFactors = FALSE)
p2_df <- read.csv(paste0(results_folder,melt_files[grep("exp2_testStatistics",melt_files)]), stringsAsFactors = FALSE)
p3_df <- read.csv(paste0(results_folder,melt_files[grep("exp3_testStatistics",melt_files)]), stringsAsFactors = FALSE)
# bootstrap dataframe contains information about the p values (obtained for tree proportion using a parametric bootstrap)
bs_df <- read.csv(paste0(results_folder,melt_files[grep("exp3_p_value",melt_files)]), stringsAsFactors = FALSE)

# Code for a simple exploratory plot - grids test statistic and tree age, each small plot is event type against test statistic value
# d = read.csv("~/Dropbox/Projects_Current/tree_likeness/results/plot1_testStatistics_collatedSimulationData_melted.csv")
# #e = subset(d, tree_age == 0.1)
# e = subset(d, tree1_tree_shape == 'balanced')
# e = subset(e, tree2_tree_shape == 'balanced')
# e$type = paste(e$tree2_event_type, e$tree2_event_position)
# ggplot(e, aes(x = type, y = value)) +
#   geom_boxplot() +
#   facet_grid(variable~tree_age, scales = "free_y") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Plots for thesis/seminar
# Plot 1: How do different events impact detection of recombination?
print("Plot 1")
e = subset(plot1_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1 )
e = subset(e, variable == "PHI_observed" | variable == "splittable_percentage"  | variable == "pdm_difference" | variable == "proportion_recombinant_triplets" |
             variable == "neighbour_net_untrimmed" | variable == "neighbour_net_trimmed" | variable == "prop_resolved_quartets" | variable == "mean_delta_q" |
             variable == "mode_delta_q")
e$type = paste(e$tree2_event_type, e$tree2_event_position)
# Have to reorder variables so the gird comes out in the right way - do this using a new column that's a factor
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","pdm_difference","proportion_recombinant_triplets","neighbour_net_untrimmed",
                                       "neighbour_net_trimmed","prop_resolved_quartets","mean_delta_q","mode_delta_q"))
facet_names <- list("PHI_observed" = "PHI \n (PhiPack)","proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)",
                    "prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", "splittable_percentage" = "Distance ratio \n (This paper)",
                    "pdm_difference" = "Distance difference \n (This paper)","neighbour_net_untrimmed" = "Tree proportion \n (Untrimmed) \n (This paper)",
                    "neighbour_net_trimmed" = "Tree proportion \n (Trimmed) \n (This paper)", "mean_delta_q" = "Mean delta_q \n (delta plots)",
                    "mode_delta_q" = "Mode delta_q \n (delta plots)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = type, y = value)) +
  geom_boxplot(outlier.size = 3) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), ncol=3) +
  scale_x_discrete(name = "\n Type of introgression event \n",
                   labels=c("none none" = "None", "reciprocal close" = "Reciprocal, \n Close", "reciprocal divergent" = "Reciprocal, \n Divergent", "reciprocal ancient" = "Reciprocal, \n Ancient",
                            "nonreciprocal close" = "Nonreciprocal, \n Close", "nonreciprocal divergent" = "Nonreciprocal, \n Divergent", "nonreciprocal ancient" = "Nonreciprocal, \n Ancient"),
                   limits=c("none none","reciprocal close","nonreciprocal close","reciprocal divergent","nonreciprocal divergent","reciprocal ancient","nonreciprocal ancient")) +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")))
# To get proper test size etc, save with the following dimensions: 4090 x 1938
ggsave(filename = paste0(plots_folder,"plot1_differentEventTypes.png"), plot = p, units = "in", width = 43, height = 43)


# Plot 2: How does increasing the proportion of the recombinant sequence affect detection of treelikeness?
print("Plot 2")
e = subset(plot2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e = subset(e, variable == "PHI_observed" | variable == "splittable_percentage"  | variable == "pdm_difference" | variable == "proportion_recombinant_triplets" |
             variable == "neighbour_net_untrimmed" | variable == "neighbour_net_trimmed" | variable == "prop_resolved_quartets" | variable == "mean_delta_q" |
             variable == "mode_delta_q")
e$type = paste(e$tree2_event_type, e$tree2_event_position)
# Have to reorder variables so the gird comes out in the right way - do this using a new column that's a factor
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","pdm_difference","proportion_recombinant_triplets","neighbour_net_untrimmed",
                                       "neighbour_net_trimmed","prop_resolved_quartets","mean_delta_q","mode_delta_q"))
facet_names <- list("PHI_observed" = "PHI \n (PhiPack)","proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)",
                    "prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", "splittable_percentage" = "Distance ratio \n (This paper)",
                    "pdm_difference" = "Distance difference \n (This paper)","neighbour_net_untrimmed" = "Tree proportion \n (Untrimmed) \n (This paper)",
                    "neighbour_net_trimmed" = "Tree proportion \n (Trimmed) \n (This paper)", "mean_delta_q" = "Mean delta_q \n (delta plots)",
                    "mode_delta_q" = "Mode delta_q \n (delta plots)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = proportion_tree2, y = value)) +
  geom_point() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), ncol=3) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")))
ggsave(filename = paste0(plots_folder,"plot2_increasingProportionTree2.png"), plot = p, units = "in", width = 43, height = 43)

# Plot 3: How does tree age affect detection of treelikeness?
print("Plot 3")
e = subset(plot2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$age = factor(e$tree_age)
e = subset(e, variable == "PHI_observed" | variable == "splittable_percentage"  | variable == "pdm_difference" | variable == "proportion_recombinant_triplets" |
             variable == "neighbour_net_untrimmed" | variable == "neighbour_net_trimmed" | variable == "prop_resolved_quartets" | variable == "mean_delta_q" |
             variable == "mode_delta_q")
# Have to reorder variables so the gird comes out in the right way - do this using a new column that's a factor
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","pdm_difference","proportion_recombinant_triplets","neighbour_net_untrimmed",
                                       "neighbour_net_trimmed","prop_resolved_quartets","mean_delta_q","mode_delta_q"))
facet_names <- list("PHI_observed" = "PHI \n (PhiPack)","proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)",
                    "prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", "splittable_percentage" = "Distance ratio \n (This paper)",
                    "pdm_difference" = "Distance difference \n (This paper)","neighbour_net_untrimmed" = "Tree proportion \n (Untrimmed) \n (This paper)",
                    "neighbour_net_trimmed" = "Tree proportion \n (Trimmed) \n (This paper)", "mean_delta_q" = "Mean delta_q \n (delta plots)",
                    "mode_delta_q" = "Mode delta_q \n (delta plots)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}

colour_list <- c("1" = "#000000", "0.5" = "#E69F00", "0.1" = "#56B4E9", "0.05" = "#CC79A7")
p <- ggplot(e, aes(x = proportion_tree2, y = value, color = age )) +
  geom_smooth(size = 3) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 3, ncol = 3) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(size = 45), axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 45), strip.text = element_text(size = 60), legend.text = element_text(size = 50),
        legend.title = element_text(size = 60), legend.key.width = unit(4,"cm"), legend.key.height = unit(2, "cm"),
        strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")), panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "#999999"),panel.grid.minor = element_line(colour = "grey78")) +
  scale_color_manual(values = colour_list) +
  guides(color = guide_legend(title = "Tree depth"))
ggsave(filename = paste0(plots_folder,"plot3_treeAgeWithIncreasingTree2.png"), plot = p, units = "in", width = 67.5, height = 46.8, limitsize = FALSE)

# Plot 4: How does the number of events impact detection of tree likeness?
print("Plot 4")
e = subset(plot3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "reciprocal")
e$event_asfactor <- as.factor(e$number_of_events)
e = subset(e, variable == "PHI_observed" | variable == "splittable_percentage"  | variable == "pdm_difference" | variable == "proportion_recombinant_triplets" |
             variable == "neighbour_net_untrimmed" | variable == "neighbour_net_trimmed" | variable == "prop_resolved_quartets" | variable == "mean_delta_q" |
             variable == "mode_delta_q")
# Have to reorder variables so the gird comes out in the right way - do this using a new column that's a factor
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","pdm_difference","proportion_recombinant_triplets","neighbour_net_untrimmed",
                                       "neighbour_net_trimmed","prop_resolved_quartets","mean_delta_q","mode_delta_q"))
facet_names <- list("PHI_observed" = "PHI \n (PhiPack)","proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)",
                    "prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", "splittable_percentage" = "Distance ratio \n (This paper)",
                    "pdm_difference" = "Distance difference \n (This paper)","neighbour_net_untrimmed" = "Tree proportion \n (Untrimmed) \n (This paper)",
                    "neighbour_net_trimmed" = "Tree proportion \n (Trimmed) \n (This paper)", "mean_delta_q" = "Mean delta_q \n (delta plots)",
                    "mode_delta_q" = "Mode delta_q \n (delta plots)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}

p <- ggplot(e, aes(x = event_asfactor, y = value)) +
  geom_boxplot(outlier.size = 5, lwd = 2) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 3, ncol = 3) +
  scale_x_discrete(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(size = 45), axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 45), strip.text = element_text(size = 60), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        panel.background = element_rect(fill="white"),panel.grid.major = element_line(colour = "grey93"),panel.grid.minor = element_line(colour = "grey98"))
ggsave(filename = paste0(plots_folder,"plot4_numberOfEvents.png"), plot = p, units = "in", width = 60, height = 46.8, limitsize = FALSE)

# Plot 5: How does reciprocity of events influence detection of treelikeness?
print("Plot 5")
e = subset(plot3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, variable == "PHI_observed" | variable == "splittable_percentage"  | variable == "pdm_difference" | variable == "proportion_recombinant_triplets" |
             variable == "neighbour_net_untrimmed" | variable == "neighbour_net_trimmed" | variable == "prop_resolved_quartets" | variable == "mean_delta_q" |
             variable == "mode_delta_q")
# Have to reorder variables so the gird comes out in the right way - do this using a new column that's a factor
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","pdm_difference","proportion_recombinant_triplets","neighbour_net_untrimmed",
                                       "neighbour_net_trimmed","prop_resolved_quartets","mean_delta_q","mode_delta_q"))
facet_names <- list("PHI_observed" = "PHI \n (PhiPack)","proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)",
                    "prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", "splittable_percentage" = "Distance ratio \n (This paper)",
                    "pdm_difference" = "Distance difference \n (This paper)","neighbour_net_untrimmed" = "Tree proportion \n (Untrimmed) \n (This paper)",
                    "neighbour_net_trimmed" = "Tree proportion \n (Trimmed) \n (This paper)", "mean_delta_q" = "Mean delta_q \n (delta plots)",
                    "mode_delta_q" = "Mode delta_q \n (delta plots)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}

p <- ggplot(e, aes(x = number_of_events, y = value, colour = tree2_event_type)) +
  geom_smooth(size = 3) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 3, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n", breaks = c(0:8), labels = c(0:8)) +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(size = 45), axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 45), strip.text = element_text(size = 60), legend.text = element_text(size = 45),
        legend.title = element_text(size = 60), legend.key.width = unit(4,"cm"), legend.key.height = unit(2, "cm"), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        panel.background = element_rect(fill="white"),panel.grid.major = element_line(colour = "#999999"),panel.grid.minor = element_line(colour = "grey78")) + 
  guides(color = guide_legend(title = "Event type")) +
  scale_colour_manual(values = c("#ca0020","#0571b0"),labels = c("Nonreciprocal", "Reciprocal"))
ggsave(filename = paste0(plots_folder,"plot5_ReciprocialAndNonreciprocalEvents.png"), plot = p, units = "in", width = 60, height = 46.8, limitsize = FALSE)

# Plot 6: Are the results statistically significant?
print("Plot 6")
e = subset(bs_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e = subset(e, variable != "PHI_observed_p_value")
e = subset(e, variable != "num_recombinant_sequences_p_value")

e = subset(e, variable == "PHI_p_value" | variable == "splittable_percentage_p_value"  | variable == "pdm_difference_p_value" | variable == "X3Seq_p_value" |
             variable == "neighbour_net_trimmed_p_value" | variable == "neighbour_net_untrimmed_p_value" | variable == "likelihood_mapping_p_value" | variable == "mean_delta_q_p_value" |
             variable == "mode_delta_q_p_value")
# Have to reorder variables so the gird comes out in the right way - do this using a new column that's a factor
e$group = factor(e$variable,levels = c("PHI_p_value","splittable_percentage_p_value","pdm_difference_p_value","X3Seq_p_value","neighbour_net_trimmed_p_value","neighbour_net_untrimmed_p_value",
                   "likelihood_mapping_p_value","mean_delta_q_p_value","mode_delta_q_p_value"))
facet_names <- list("PHI_p_value" = "PHI \n (PhiPack)","splittable_percentage_p_value" = "Distance \n ratio","pdm_difference_p_value" = "Distance \n difference",
                   "X3Seq_p_value" = "Proportion of \n recombinant triplets \n (3SEQ)","neighbour_net_trimmed_p_value" = "Tree \n proportion \n (Trimmed)",
                   "neighbour_net_untrimmed_p_value" = "Tree \n proportion \n (Untrimmed)","likelihood_mapping_p_value" = "Proportion of \n resolved quartets \n (IQ-Tree)",
                   "mean_delta_q_p_value" = "Mean delta q","mode_delta_q_p_value" = "Mode delta q")

facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~proportion_tree2,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  xlab("\n P value \n") +
  ylab("\n Count \n") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        strip.text.y = element_text(margin = margin(1,0,0.5,0, "cm"))) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1))
ggsave(filename = paste0(plots_folder,"plot6_StatisticalSignificance.png"), plot = p, units = "in", width = 43, height = 63, limitsize = FALSE)

# Plot 7: Are the tree proportion results statistically significant?
print("Plot 7")
e = subset(bs_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e = subset(e, variable != "PHI_observed_p_value")
e = subset(e, variable != "num_recombinant_sequences_p_value")
e = subset(e, variable == "neighbour_net_trimmed_p_value" | variable == "neighbour_net_untrimmed_p_value")
e$group = factor(e$variable,levels = c("neighbour_net_untrimmed_p_value","neighbour_net_trimmed_p_value"))
facet_names <- list("neighbour_net_untrimmed_p_value" = "Tree proportion (Untrimmed)","neighbour_net_trimmed_p_value" = "Tree proportion (Trimmed)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~proportion_tree2,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  xlab("\n P value \n") +
  ylab("\n Count \n") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        strip.text.y = element_text(margin = margin(1,0,0.5,0, "cm"))) +
        scale_x_continuous(labels = c(0,0.25,0.5,0.75,1))
ggsave(filename = paste0(plots_folder,"plot7_StatisticalSignificance_treeProportion.png"), plot = p, units = "in", width = 43, height = 22)

# Plot 8
print("Plot 8")
e = subset(bs_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e = subset(e, variable != "PHI_observed_p_value")
e = subset(e, variable != "num_recombinant_sequences_p_value")
# to make new df for the plot
PHI_p_value <- c(nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                 nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
X3SEQ_p_value <- c(nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                   nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
likelihood_mapping_p_value <- c(nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                                nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
splittable_percentage_p_value <- c(nrow(e[e$variable == "splittable_percentage_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "splittable_percentage_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "splittable_percentage_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                                   nrow(e[e$variable == "splittable_percentage_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "splittable_percentage_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "splittable_percentage_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
pdm_difference_p_value <- c(nrow(e[e$variable == "pdm_difference_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "pdm_difference_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "pdm_difference_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                             nrow(e[e$variable == "pdm_difference_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "pdm_difference_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "pdm_difference_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
neighbour_net_trimmed_p_value <- c(nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                                  nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
neighbour_net_untrimmed_p_value <- c(nrow(e[e$variable == "neighbour_net_untrimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "neighbour_net_untrimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "neighbour_net_untrimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                                     nrow(e[e$variable == "neighbour_net_untrimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "neighbour_net_untrimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "neighbour_net_untrimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
mean_delta_q_p_value <- c(nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                          nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
mode_delta_q_p_value <- c(nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                          nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))


value <- c(PHI_p_value,X3SEQ_p_value,likelihood_mapping_p_value,splittable_percentage_p_value,pdm_difference_p_value,neighbour_net_untrimmed_p_value,neighbour_net_trimmed_p_value,mean_delta_q_p_value,mode_delta_q_p_value)
ts <- c(rep("PHI_p_value",6), rep("X3SEQ_p_value",6), rep("likelihood_mapping_p_value",6), rep("splittable_percentage_p_value",6), rep("pdm_difference_p_value",6), rep("neighbour_net_untrimmed_p_value",6), rep("neighbour_net_trimmed_p_value",6), rep("mean_delta_q_p_value",6), rep("mode_delta_q_p_value",6))
proportion_introgressed_DNA <- c(rep(seq(0,0.5,0.1),9))
f <- data.frame(proportion_introgressed_DNA,ts,value,stringsAsFactors = FALSE)
f$variable <- factor(ts, levels = c("PHI_p_value","splittable_percentage_p_value","pdm_difference_p_value","X3SEQ_p_value","neighbour_net_untrimmed_p_value",
                                    "neighbour_net_trimmed_p_value","likelihood_mapping_p_value","mean_delta_q_p_value","mode_delta_q_p_value"))
facet_names <- list("PHI_p_value" = "PHI \n (PhiPack)","splittable_percentage_p_value" = "Distance \n ratio","pdm_difference_p_value" = "Distance \n difference",
                    "X3SEQ_p_value" = "Proportion of \n recombinant triplets \n (3SEQ)","neighbour_net_trimmed_p_value" = "Tree \n proportion \n (Trimmed)",
                    "neighbour_net_untrimmed_p_value" = "Tree \n proportion \n (Untrimmed)","likelihood_mapping_p_value" = "Proportion of \n resolved quartets \n (IQ-Tree)",
                    "mean_delta_q_p_value" = "Mean delta q","mode_delta_q_p_value" = "Mode delta q")

facet_labeller <- function(variable){
  variable <- facet_names[variable]
}

p <- ggplot(f, aes(x = proportion_introgressed_DNA, y = value)) +
  geom_line(size=3) +
  facet_wrap(~variable, labeller = labeller(variable = facet_labeller), nrow = 3, ncol = 3) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Percent of simulations that reject the null hypothesis \n (p-value < 0.05) \n") +
  theme(axis.text.x = element_text(size = 40), axis.title.x = element_text(size = 50), axis.title.y = element_text(size = 50),
        axis.text.y = element_text(size = 40), strip.text = element_text(size = 50), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        strip.text.y = element_text(margin = margin(1,0,0.5,0, "cm")), legend.text = element_text(size = 40),
        legend.title = element_text(size = 50), legend.key.width = unit(5,"cm"), legend.key.height = unit(5, "cm"),
        panel.background = element_rect(fill="white"),panel.grid.major = element_line(colour = "#999999"),panel.grid.minor = element_line(colour = "grey78")) +
  scale_y_continuous(labels = seq(0,100,10), breaks = seq(0,100,10), minor_breaks = seq(0,100,5), limits = c(0,100)) + 
  geom_hline(aes(yintercept = 5, colour = "red"), linetype = "dashed", size = 1.5) + 
  scale_colour_manual("Ideal false\npositive rate\n", values="red", labels = "5% when\n\u03b1 = 0.05")
ggsave(filename = paste0(plots_folder,"plot8_facetedPValues_fixedAxes.png"), plot = p, units = "in", width = 40, height = 46.8)
