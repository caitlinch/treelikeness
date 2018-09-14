# R code to plot simulation results
# Specify which file paths to use

run_location = "mac"
# run_location = "soma"

if (run_location == "mac"){
  # Set file paths etc
  input_folder <- "/Users/caitlincherryh/Documents/Results/collatedOutput/"
  output_folder <- "/Users/caitlincherryh/Documents/Results/plots/seminarPlots/attempt_2/"
  # Set working directory
  maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # for work computer
} else if (run_location == "soma") {
  # Set file paths etc
  input_folder <- "/data/caitlin/treelikeness/results/"
  
  # Set working directory
  maindir <- "/data/caitlin/treelikeness/"
}

# load required libraries
library(ggplot2)
library(gridExtra)
library(viridis)

# Open dataframes
plot1_df <- read.csv(paste0(input_folder,"plot1_testStatistics_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)
plot2_df <- read.csv(paste0(input_folder,"plot2_testStatistics_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)
plot3_df <- read.csv(paste0(input_folder,"plot3_testStatistics_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)
plot4_df <- read.csv(paste0(input_folder,"plot4_testStatistics_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)
bs_df <-  read.csv(paste0(input_folder,"plot4_p_value_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)

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
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree-splits proportion \n (Neighbor-Net)","pdm_difference" = "Total element-wise difference","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Modified splittable percentage")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = type, y = value)) +
  geom_boxplot() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_discrete(name = "Type of introgression event",
                   labels=c("none none" = "None", "reciprocal close" = "Reciprocal, \n Close", "reciprocal divergent" = "Reciprocal, \n Divergent", "reciprocal ancient" = "Reciprocal, \n Ancient",
                            "nonreciprocal close" = "Nonreciprocal, \n Close", "nonreciprocal divergent" = "Nonreciprocal, \n Divergent", "nonreciprocal ancient" = "Nonreciprocal, \n Ancient"),
                   limits=c("none none","reciprocal close","nonreciprocal close","reciprocal divergent","nonreciprocal divergent","reciprocal ancient","nonreciprocal ancient")) +
  ylab("Test statistic value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")))
# To get proper test size etc, save with the following dimensions: 4090 x 1938
ggsave(filename = paste0(output_folder,"plot1_differentEventTypes.png"), plot = p, units = "in", width = 43, height = 20.4)


# Plot 2: How does increasing the proportion of the recombinant sequence affect detection of treelikeness?
print("Plot 2")
e = subset(plot2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree-splits proportion \n (Neighbor-Net)","pdm_difference" = "Total element-wise difference","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Modified splittable percentage")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = proportion_tree2, y = value)) +
  geom_point() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "Proportion of tree 2") +
  ylab("Test statistic value") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")))
ggsave(filename = paste0(output_folder,"plot2_increasingProportionTree2.png"), plot = p, units = "in", width = 43, height = 20.4)

# Plot 3: How does tree age affect detection of treelikeness?
print("Plot 3")
e = subset(plot2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$age = factor(e$tree_age)
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree-splits proportion \n (Neighbor-Net)","pdm_difference" = "Total element-wise difference","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Modified splittable percentage")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = proportion_tree2, y = value, color = age )) +
  geom_smooth(size = 3) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "Proportion of tree 2") +
  ylab("Test statistic value") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), legend.text = element_text(size = 30),
        legend.title = element_text(size = 40), legend.key.width = unit(4,"cm"), legend.key.height = unit(2, "cm"),
        strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm"))) +
  scale_color_viridis(discrete = TRUE, option = "viridis") +
  guides(color = guide_legend(title = "Tree depth"))
ggsave(filename = paste0(output_folder,"plot3_treeAgeWithIncreasingTree2.png"), plot = p, units = "in", width = 48, height = 20.4)

# Plot 4: How does the number of events impact detection of tree likeness?
print("Plot 4")
e = subset(plot3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$event_asfactor <- as.factor(e$number_of_events)
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree-splits proportion \n (Neighbor-Net)","pdm_difference" = "Total element-wise difference","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Modified splittable percentage")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = event_asfactor, y = value)) +
  geom_boxplot() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_discrete(name = "Number of introgression events") +
  ylab("Test statistic value") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")))
ggsave(filename = paste0(output_folder,"plot4_numberOfEvents.png"), plot = p, units = "in", width = 43, height = 20.4)

# Plot 5: How does reciprocity of events influence detection of treelikeness?
print("Plot 5")
e = subset(plot3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree-splits proportion \n (Neighbor-Net)","pdm_difference" = "Total element-wise difference","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Modified splittable percentage")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = number_of_events, y = value, colour = tree2_event_type)) +
  geom_smooth(size = 3) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "Number of introgression events", breaks = c(0:8), labels = c(0:8)) +
  ylab("Test statistic value") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), legend.text = element_text(size = 30),
        legend.title = element_text(size = 40), legend.key.width = unit(4,"cm"), legend.key.height = unit(2, "cm"), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm"))) + 
  scale_color_viridis(discrete = TRUE, option = "plasma", labels = c("Nonreciprocal", "Reciprocal")) + 
  guides(color = guide_legend(title = "Event type"))
ggsave(filename = paste0(output_folder,"plot5_ReciprocialAndNonreciprocalEvents.png"), plot = p, units = "in", width = 48, height = 20.4)

# Plot 6: Are the results statistically significant?
print("Plot 6")
e = subset(bs_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$group = factor(e$variable,levels = c("PHI_p_value","X3Seq_p_value","likelihood_mapping_p_value","splittable_percentage_p_value","pdm_difference_p_value","neighbour_net_p_value"))
facet_names <- list("neighbour_net_p_value" = "Tree-splits \n proportion \n (Neighbor-Net)","pdm_difference_p_value" = "Total element-wise \n difference","PHI_p_value" = "PHI \n (PhiPack)","likelihood_mapping_p_value" = "Proportion of \n resolved quartets \n (IQ-Tree)", 
                    "X3Seq_p_value" = "Proportion of \n recombinant triplets \n (3SEQ)", "splittable_percentage_p_value" = "Modified splittable \n percentage")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~proportion_tree2,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "P value") +
  ylab("Count") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        strip.text.y = element_text(margin = margin(1,0,0.5,0, "cm"))) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1))
ggsave(filename = paste0(output_folder,"plot6_StatisticalSignificance.png"), plot = p, units = "in", width = 43, height = 43)

