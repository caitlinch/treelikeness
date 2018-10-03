# R code to plot simulation results
# Specify which file paths to use

run_location = "mac"
# run_location = "soma"

if (run_location == "mac"){
  # Set file paths etc
  input_folder <- "/Users/caitlincherryh/Documents/Results/simulations_20180913/collatedOutput/"
  output_folder <- "/Users/caitlincherryh/Documents/Results/simulations_20180913/plots/plots_20181003/"
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
facet_names <- list("neighbour_net" = "Tree proportion","pdm_difference" = "Distance difference","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Distance ratio")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = type, y = value)) +
  geom_boxplot() +
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


# Plot 2: How does increasing the proportion of the recombinant sequence affect detection of treelikeness?
print("Plot 2")
e = subset(plot2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree proportion","pdm_difference" = "Distance difference","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Distance ratio")
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

# Plot 3: How does tree age affect detection of treelikeness?
print("Plot 3")
e = subset(plot2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$age = factor(e$tree_age)
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree proportion","pdm_difference" = "Distance difference","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Distance ratio")
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

# Plot 4: How does the number of events impact detection of tree likeness?
print("Plot 4")
e = subset(plot3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "reciprocal")
e$event_asfactor <- as.factor(e$number_of_events)
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree proportion","pdm_difference" = "Distance difference","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Distance ratio")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = event_asfactor, y = value)) +
  geom_boxplot() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_discrete(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")))
ggsave(filename = paste0(output_folder,"plot4_numberOfEvents.png"), plot = p, units = "in", width = 43, height = 20.4)

# Plot 5: How does reciprocity of events influence detection of treelikeness?
print("Plot 5")
e = subset(plot3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree proportion","pdm_difference" = "Distance difference","PHI_observed" = "PHI \n (PhiPack)","prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)", "splittable_percentage" = "Distance ratio")
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

# Plot 6: Are the results statistically significant?
print("Plot 6")
e = subset(bs_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
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

