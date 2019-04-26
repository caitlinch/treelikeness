# R code to plot simulation results
# Specify which file paths to use

run_location = "mac"
# run_location = "soma"

if (run_location == "mac"){
  # Set file paths etc
  input_folder <- "/Users/caitlincherryh/Documents/Honours/Results/simulations_20190411/collatedOutput/"
  output_folder <- "/Users/caitlincherryh/Documents/Honours/Results/simulations_20190411/plots/"
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
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","proportion_recombinant_triplets","pdm_difference","prop_resolved_quartets","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree proportion \n (This thesis)","pdm_difference" = "Distance difference \n (This thesis)","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Distance ratio \n (This thesis)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = type, y = value)) +
  geom_boxplot(outlier.size = 3) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_discrete(name = "\n Type of introgression event \n",
                   labels=c("none none" = "None", "reciprocal close" = "Reciprocal, \n Close", "reciprocal divergent" = "Reciprocal, \n Divergent", "reciprocal ancient" = "Reciprocal, \n Ancient",
                            "nonreciprocal close" = "Nonreciprocal, \n Close", "nonreciprocal divergent" = "Nonreciprocal, \n Divergent", "nonreciprocal ancient" = "Nonreciprocal, \n Ancient"),
                   limits=c("none none","reciprocal close","nonreciprocal close","reciprocal divergent","nonreciprocal divergent","reciprocal ancient","nonreciprocal ancient")) +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")))
# To get proper test size etc, save with the following dimensions: 4090 x 1938
ggsave(filename = paste0(output_folder,"plot1_differentEventTypes.png"), plot = p, units = "in", width = 43, height = 20.4)

p <- ggplot(e, aes(x = type, y = value)) +
  geom_boxplot(outlier.size = 5, lwd = 2) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 3, ncol = 2) +
  scale_x_discrete(name = "\n Type of introgression event \n",
                   labels=c("none none" = "None", "reciprocal close" = "Reciprocal, \n Close", "reciprocal divergent" = "Reciprocal, \n Divergent", "reciprocal ancient" = "Reciprocal, \n Ancient",
                            "nonreciprocal close" = "Nonreciprocal, \n Close", "nonreciprocal divergent" = "Nonreciprocal, \n Divergent", "nonreciprocal ancient" = "Nonreciprocal, \n Ancient"),
                   limits=c("none none","reciprocal close","nonreciprocal close","reciprocal divergent","nonreciprocal divergent","reciprocal ancient","nonreciprocal ancient")) +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 45), axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 45), strip.text = element_text(size = 60), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        panel.background = element_rect(fill="white"),panel.grid.major = element_line(colour = "grey93"),panel.grid.minor = element_line(colour = "grey98"))
# To get proper test size etc, save with the following dimensions: 4090 x 1938
ggsave(filename = paste0(output_folder,"plot1_differentEventTypes_portrait.png"), plot = p, units = "in", width = 40, height = 46.8)


# Plot 2: How does increasing the proportion of the recombinant sequence affect detection of treelikeness?
print("Plot 2")
e = subset(plot2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","proportion_recombinant_triplets","pdm_difference","prop_resolved_quartets","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree proportion \n (This thesis)","pdm_difference" = "Distance difference \n (This thesis)","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Distance ratio \n (This thesis)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = proportion_tree2, y = value)) +
  geom_point() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")))
ggsave(filename = paste0(output_folder,"plot2_increasingProportionTree2.png"), plot = p, units = "in", width = 43, height = 20.4)

p <- ggplot(e, aes(x = proportion_tree2, y = value)) +
  geom_point(size=5) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 3, ncol = 2) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(size = 45), axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 45), strip.text = element_text(size = 60), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        panel.background = element_rect(fill="white"),panel.grid.major = element_line(colour = "#999999"),panel.grid.minor = element_line(colour = "grey78"))
ggsave(filename = paste0(output_folder,"plot2_increasingProportionTree2_portrait.png"), plot = p, units = "in", width = 40, height = 46.8)

# Plot 3: How does tree age affect detection of treelikeness?
print("Plot 3")
e = subset(plot2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$age = factor(e$tree_age)
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","proportion_recombinant_triplets","pdm_difference","prop_resolved_quartets","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree proportion \n (This thesis)","pdm_difference" = "Distance difference \n (This thesis)","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Distance ratio \n (This thesis)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = proportion_tree2, y = value, color = age )) +
  geom_smooth(size = 3) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), legend.text = element_text(size = 30),
        legend.title = element_text(size = 40), legend.key.width = unit(4,"cm"), legend.key.height = unit(2, "cm"),
        strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm"))) +
  scale_color_viridis(discrete = TRUE, option = "viridis") +
  guides(color = guide_legend(title = "Tree depth"))
ggsave(filename = paste0(output_folder,"plot3_treeAgeWithIncreasingTree2.png"), plot = p, units = "in", width = 48, height = 20.4)

colour_list <- c("1" = "#000000", "0.5" = "#E69F00", "0.1" = "#56B4E9", "0.05" = "#CC79A7")
p <- ggplot(e, aes(x = proportion_tree2, y = value, color = age )) +
  geom_smooth(size = 3) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 3, ncol = 2) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(size = 45), axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 45), strip.text = element_text(size = 60), legend.text = element_text(size = 50),
        legend.title = element_text(size = 60), legend.key.width = unit(4,"cm"), legend.key.height = unit(2, "cm"),
        strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")), panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "#999999"),panel.grid.minor = element_line(colour = "grey78")) +
  scale_color_manual(values = colour_list) +
  guides(color = guide_legend(title = "Tree depth"))
ggsave(filename = paste0(output_folder,"plot3_treeAgeWithIncreasingTree2_portrait.png"), plot = p, units = "in", width = 45, height = 46.8)

# Plot 4: How does the number of events impact detection of tree likeness?
print("Plot 4")
e = subset(plot3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "reciprocal")
e$event_asfactor <- as.factor(e$number_of_events)
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","proportion_recombinant_triplets","pdm_difference","prop_resolved_quartets","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree proportion \n (This thesis)","pdm_difference" = "Distance difference \n (This thesis)","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Distance ratio \n (This thesis)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = event_asfactor, y = value)) +
  geom_boxplot(outlier.size = 3) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_discrete(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")))
ggsave(filename = paste0(output_folder,"plot4_numberOfEvents.png"), plot = p, units = "in", width = 43, height = 20.4)

p <- ggplot(e, aes(x = event_asfactor, y = value)) +
  geom_boxplot(outlier.size = 5, lwd = 2) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 3, ncol = 2) +
  scale_x_discrete(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(size = 45), axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 45), strip.text = element_text(size = 60), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        panel.background = element_rect(fill="white"),panel.grid.major = element_line(colour = "grey93"),panel.grid.minor = element_line(colour = "grey98"))
ggsave(filename = paste0(output_folder,"plot4_numberOfEvents_portrait.png"), plot = p, units = "in", width = 40, height = 46.8)

# Plot 5: How does reciprocity of events influence detection of treelikeness?
print("Plot 5")
e = subset(plot3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","proportion_recombinant_triplets","pdm_difference","prop_resolved_quartets","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree proportion \n (This thesis)","pdm_difference" = "Distance difference \n (This thesis)","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Distance ratio \n (This thesis)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = number_of_events, y = value, colour = tree2_event_type)) +
  geom_smooth(size = 3) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "\n Number of introgression events \n", breaks = c(0:8), labels = c(0:8)) +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), legend.text = element_text(size = 30),
        legend.title = element_text(size = 40), legend.key.width = unit(4,"cm"), legend.key.height = unit(2, "cm"), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm"))) + 
  scale_color_viridis(discrete = TRUE, option = "plasma", labels = c("Nonreciprocal", "Reciprocal")) + 
  guides(color = guide_legend(title = "Event type"))
ggsave(filename = paste0(output_folder,"plot5_ReciprocialAndNonreciprocalEvents.png"), plot = p, units = "in", width = 48, height = 20.4)

p <- ggplot(e, aes(x = number_of_events, y = value, colour = tree2_event_type)) +
  geom_smooth(size = 3) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 3, ncol = 2) +
  scale_x_continuous(name = "\n Number of introgression events \n", breaks = c(0:8), labels = c(0:8)) +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(size = 45), axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 45), strip.text = element_text(size = 60), legend.text = element_text(size = 45),
        legend.title = element_text(size = 60), legend.key.width = unit(4,"cm"), legend.key.height = unit(2, "cm"), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        panel.background = element_rect(fill="white"),panel.grid.major = element_line(colour = "#999999"),panel.grid.minor = element_line(colour = "grey78")) + 
  guides(color = guide_legend(title = "Event type")) +
  scale_colour_manual(values = c("#ca0020","#0571b0"),labels = c("Nonreciprocal", "Reciprocal"))
ggsave(filename = paste0(output_folder,"plot5_ReciprocialAndNonreciprocalEvents_portrait.png"), plot = p, units = "in", width = 40, height = 46.8)

# Plot 6: Are the results statistically significant?
print("Plot 6")
e = subset(bs_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e = subset(e, variable != "PHI_observed_p_value")
e = subset(e, variable != "num_recombinant_sequences_p_value")
e$group = factor(e$variable,levels = c("PHI_p_value","X3Seq_p_value","likelihood_mapping_p_value","splittable_percentage_p_value","pdm_difference_p_value","neighbour_net_p_value"))
facet_names <- list("neighbour_net_p_value" = "Tree \n proportion","pdm_difference_p_value" = "Distance \n difference","PHI_p_value" = "PHI \n (PhiPack)","likelihood_mapping_p_value" = "Proportion of \n resolved quartets \n (IQ-Tree)", 
                    "X3Seq_p_value" = "Proportion of \n recombinant triplets \n (3SEQ)", "splittable_percentage_p_value" = "Distance \n ratio")
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
ggsave(filename = paste0(output_folder,"plot6_StatisticalSignificance.png"), plot = p, units = "in", width = 43, height = 43)

# Plot 7: Are the tree proportion results statistically significant?
print("Plot 7")
e = subset(bs_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e = subset(e, variable != "PHI_observed_p_value")
e = subset(e, variable != "num_recombinant_sequences_p_value")
e = subset(e, variable == "neighbour_net_p_value")
e$group = factor(e$variable,levels = c("neighbour_net_p_value"))
facet_names <- list("neighbour_net_p_value" = "Tree proportion")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(~proportion_tree2,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  xlab("\n P value \n") +
  ylab("\n Count \n") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        strip.text.y = element_text(margin = margin(1,0,0.5,0, "cm"))) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1))
ggsave(filename = paste0(output_folder,"plot7_StatisticalSignificance_treeProportion.png"), plot = p, units = "in", width = 43, height = 11)

# Plot 8: Are the tree proportion results statistically significant?
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
neighbour_net_p_value <- c(nrow(e[e$variable == "neighbour_net_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "neighbour_net_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "neighbour_net_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                           nrow(e[e$variable == "neighbour_net_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "neighbour_net_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "neighbour_net_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
value <- c(PHI_p_value,X3SEQ_p_value,likelihood_mapping_p_value,splittable_percentage_p_value,pdm_difference_p_value,neighbour_net_p_value)
value <- c(PHI_p_value,X3SEQ_p_value,likelihood_mapping_p_value,splittable_percentage_p_value,pdm_difference_p_value,neighbour_net_p_value)
ts <- c(rep("PHI_p_value",6), rep("X3SEQ_p_value",6), rep("likelihood_mapping_p_value",6), rep("splittable_percentage_p_value",6), rep("pdm_difference_p_value",6), rep("neighbour_net_p_value",6))
proportion_introgressed_DNA <- c(rep(seq(0,0.5,0.1),6))
f <- data.frame(proportion_introgressed_DNA,ts,value,stringsAsFactors = FALSE)
f$variable <- factor(ts, levels = c("PHI_p_value","X3SEQ_p_value","likelihood_mapping_p_value","splittable_percentage_p_value","pdm_difference_p_value","neighbour_net_p_value"))

cols <- c("neighbour_net_p_value" = "#56B4E9","pdm_difference_p_value" = "#E69F00","PHI_p_value" = "#000000","likelihood_mapping_p_value" = "#009E73", 
          "X3SEQ_p_value" = "#0072B2", "splittable_percentage_p_value" = "#F0E442")
shapes <- c("neighbour_net_p_value" = 15,"pdm_difference_p_value" = 17,"PHI_p_value" = 1,"likelihood_mapping_p_value" = 19, 
            "X3SEQ_p_value" = 8, "splittable_percentage_p_value" = 18)
names <- c("neighbour_net_p_value" = "Tree \n proportion","pdm_difference_p_value" = "Distance \n difference","PHI_p_value" = "PHI \n (PhiPack)","likelihood_mapping_p_value" = "Proportion of \n resolved quartets \n (IQ-Tree)", 
          "X3SEQ_p_value" = "Proportion of \n recombinant triplets \n (3SEQ)", "splittable_percentage_p_value" = "Distance \n ratio")
p <- ggplot(f, aes(x = proportion_introgressed_DNA, y = value, color = variable, shape = variable)) +
  geom_point(size = 15) +
  geom_line(size = 3) + 
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  scale_y_continuous(name = "\n Percent of simulations that reject the null hypothesis \n (p-value < 0.05) \n", breaks = c(seq(0,100,10)), labels = c(seq(0,100,10))) +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), legend.text = element_text(size = 30),
        legend.title = element_text(size = 40), legend.key.width = unit(4,"cm"), legend.key.height = unit(4, "cm"),
        strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm"))) +
  scale_colour_manual("Test statistic", values=cols, labels = names) +
  scale_shape_manual("Test statistic", values=shapes, labels = names)
ggsave(filename = paste0(output_folder,"plot8_StatisticallySignificantResults.png"), plot = p, units = "in", width = 43, height = 20.4)

# Plot 9: facetted significance
print("Plot 9")
f$group = factor(f$variable,levels = c("PHI_p_value","splittable_percentage_p_value","X3SEQ_p_value","pdm_difference_p_value","likelihood_mapping_p_value","neighbour_net_p_value"))
facet_names <- list("neighbour_net_p_value" = "Tree proportion \n (This thesis)","pdm_difference_p_value" = "Distance difference \n (This thesis)","PHI_p_value" = "PHI \n (PhiPack)","likelihood_mapping_p_value" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "X3SEQ_p_value" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage_p_value" = "Distance ratio \n (This thesis)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(f, aes(x = proportion_introgressed_DNA, y = value)) +
  geom_line(size=3) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 3, ncol = 2) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Percent of simulations that reject the null hypothesis \n (p-value < 0.05) \n") +
  theme(axis.text.x = element_text(size = 50), axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 50), strip.text = element_text(size = 60), legend.text = element_text(size = 50),
        legend.title = element_text(size = 60), legend.key.width = unit(4,"cm"), legend.key.height = unit(2, "cm"), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        panel.background = element_rect(fill="white"),panel.grid.major = element_line(colour = "#999999"),panel.grid.minor = element_line(colour = "grey78"))
ggsave(filename = paste0(output_folder,"plot9_facetedPValues.png"), plot = p, units = "in", width = 40, height = 46.8)

# Plot 10 : parametric bootstrap p values for PHI and 3SEQ
print("Plot 10")
e = subset(bs_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
PHI_p_value <- c(nrow(e[e$variable == "PHI_observed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "PHI_observed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "PHI_observed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                 nrow(e[e$variable == "PHI_observed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "PHI_observed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "PHI_observed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
X3SEQ_p_value <- c(nrow(e[e$variable == "num_recombinant_sequences_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "num_recombinant_sequences_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "num_recombinant_sequences_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                   nrow(e[e$variable == "num_recombinant_sequences_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "num_recombinant_sequences_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "num_recombinant_sequences_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
value <- c(PHI_p_value,X3SEQ_p_value,likelihood_mapping_p_value,splittable_percentage_p_value,pdm_difference_p_value,neighbour_net_p_value)
value <- c(PHI_p_value,X3SEQ_p_value,likelihood_mapping_p_value,splittable_percentage_p_value,pdm_difference_p_value,neighbour_net_p_value)
ts <- c(rep("PHI_p_value",6), rep("X3SEQ_p_value",6), rep("likelihood_mapping_p_value",6), rep("splittable_percentage_p_value",6), rep("pdm_difference_p_value",6), rep("neighbour_net_p_value",6))
proportion_introgressed_DNA <- c(rep(seq(0,0.5,0.1),6))
f <- data.frame(proportion_introgressed_DNA,ts,value,stringsAsFactors = FALSE)
f$variable <- factor(ts, levels = c("PHI_p_value","X3SEQ_p_value","likelihood_mapping_p_value","splittable_percentage_p_value","pdm_difference_p_value","neighbour_net_p_value"))
f$group = factor(f$variable,levels = c("PHI_p_value","splittable_percentage_p_value","X3SEQ_p_value","pdm_difference_p_value","likelihood_mapping_p_value","neighbour_net_p_value"))
facet_names <- list("neighbour_net_p_value" = "Tree proportion \n (This thesis)","pdm_difference_p_value" = "Distance difference \n (This thesis)","PHI_p_value" = "PHI \n (PhiPack)","likelihood_mapping_p_value" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "X3SEQ_p_value" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage_p_value" = "Distance ratio \n (This thesis)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(f, aes(x = proportion_introgressed_DNA, y = value)) +
  geom_line(size=3) +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller), nrow = 3, ncol = 2) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  theme(axis.text.x = element_text(size = 40), axis.title.x = element_text(size = 50), axis.title.y = element_text(size = 50),
        axis.text.y = element_text(size = 40), strip.text = element_text(size = 50), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        strip.text.y = element_text(margin = margin(1,0,0.5,0, "cm")), legend.text = element_text(size = 40),
        legend.title = element_text(size = 50), legend.key.width = unit(5,"cm"), legend.key.height = unit(5, "cm"),
        panel.background = element_rect(fill="white"),panel.grid.major = element_line(colour = "#999999"),panel.grid.minor = element_line(colour = "grey78")) + 
  geom_hline(aes(yintercept = 5, colour = "red"), linetype = "dashed", size = 1.5) + 
  scale_colour_manual("Ideal false\npositive rate\n", values="red", labels = "5% when\n\u03b1 = 0.05")
ggsave(filename = paste0(output_folder,"plot10_facetedPValues_allParametric_freeAxes.png"), plot = p, units = "in", width = 40, height = 46.8)

p <- ggplot(f, aes(x = proportion_introgressed_DNA, y = value)) +
  geom_line(size=3) +
  facet_wrap(~group, labeller = labeller(group = facet_labeller), nrow = 3, ncol = 2) +
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
ggsave(filename = paste0(output_folder,"plot10_facetedPValues_allParametric_fixedAxes.png"), plot = p, units = "in", width = 40, height = 46.8)


# Plot 11 - grid of p values, all parametric
print("Plot 11")
e = subset(bs_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e = subset(e, variable != "PHI_p_value")
e = subset(e, variable != "X3Seq_p_value")
e$group = factor(e$variable,levels = c("PHI_observed_p_value","num_recombinant_sequences_p_value","likelihood_mapping_p_value","splittable_percentage_p_value","pdm_difference_p_value","neighbour_net_p_value"))
facet_names <- list("neighbour_net_p_value" = "Tree \n proportion \n (This thesis)","pdm_difference_p_value" = "Distance \n difference \n (This thesis)","PHI_observed_p_value" = "PHI \n (PhiPack)","likelihood_mapping_p_value" = "Proportion of \n resolved quartets \n (IQ-Tree)", 
                    "num_recombinant_sequences_p_value" = "Proportion of \n recombinant triplets \n (3SEQ)", "splittable_percentage_p_value" = "Distance \n ratio \n (This thesis)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}

p <- ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~proportion_tree2,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 50), axis.title.y = element_text(size = 50),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 50), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        strip.text.y = element_text(margin = margin(1,0,0.5,0, "cm")), legend.text = element_text(size = 40),
        legend.title = element_text(size = 50), legend.key.width = unit(5,"cm"), legend.key.height = unit(5, "cm"),
        panel.background = element_rect(fill="white"),panel.grid.major = element_line(colour = "#999999"),panel.grid.minor = element_line(colour = "grey78")) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1)) +
  geom_vline(aes(xintercept = 0.05, colour = "red"), linetype = "dashed", size = 1.5, show.legend = TRUE) + 
  scale_colour_manual("Statistical\nsignificance\nthreshold", values="red", labels = "p = 0.05")
ggsave(filename = paste0(output_folder,"plot11_StatisticalSignificance_allParametric_freeAxes.png"), plot = p, units = "in", width = 43, height = 43)

p <- ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~proportion_tree2, labeller = labeller(group = facet_labeller)) +
  xlab("\n P value \n") +
  ylab("\n Count \n") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 50), axis.title.y = element_text(size = 50),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 50), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        strip.text.y = element_text(margin = margin(1,0,0.5,0, "cm")), legend.text = element_text(size = 40),
        legend.title = element_text(size = 50), legend.key.width = unit(5,"cm"), legend.key.height = unit(5, "cm"),
        panel.background = element_rect(fill="white"),panel.grid.major = element_line(colour = "#999999"),panel.grid.minor = element_line(colour = "grey78")) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1)) +
  scale_y_continuous(labels = seq(0,100,20), breaks = seq(0,100,20), minor_breaks = seq(0,100,10), limits = c(0,100)) + 
  geom_vline(aes(xintercept = 0.05, colour = "red"), linetype = "dashed", size = 1.5, show.legend = TRUE) + 
  scale_colour_manual("Statistical\nsignificance\nthreshold", values="red", labels = "p = 0.05")
ggsave(filename = paste0(output_folder,"plot11_StatisticalSignificance_allParametric_fixedAxes.png"), plot = p, units = "in", width = 45, height = 43)

