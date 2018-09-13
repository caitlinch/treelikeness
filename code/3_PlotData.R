# R code to plot simulation results
# Specify which file paths to use

run_location = "mac"
# run_location = "soma"

if (run_location == "mac"){
  # Set file paths etc
  input_folder <- "/Users/caitlincherryh/Documents/Results/collatedOutput/"
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

## Exploratory plots
# Plot 1:
ggplot(data=plot1_df,aes(x=tree2_event_position,y=value)) + geom_boxplot(data=plot1_df,aes(fill = tree2_event_type)) + facet_wrap(~variable) + ylim(0,2)

# Rob plot 1
#e = subset(d, tree_age == 0.1) # = df
e = subset(plot1_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e$type = paste(e$tree2_event_type, e$tree2_event_position)
ggplot(e, aes(x = type, y = value)) +
  geom_boxplot() +
  facet_grid(variable~tree_age, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Plot 2
ggplot(data=plot2_df,aes(x=proportion_tree2,y=value)) + geom_point(data=plot2_df, aes(colour = as.factor(tree_age))) + facet_wrap(~variable) + geom_smooth(method = "lm", data = plot2_df, aes(x=proportion_tree2,y=value, colour = as.factor(tree_age)))
ggplot(data=plot2_df,aes(x=proportion_tree2,y=value)) + facet_wrap(~variable) + geom_smooth(method = "lm", data = plot2_df, aes(x=proportion_tree2,y=value, colour = as.factor(tree_age)))
# Plot 2 attempt 2
# reciprocal balanced
df <- subset(plot2_df, tree1_tree_shape == 'balanced')
df <- subset(df, tree2_tree_shape == 'balanced')
df <- subset(df, tree2_event_type == 'reciprocal')
df <- subset(df, tree2_event_position == 'close')
df$type = paste(df$tree2_event_type, df$tree2_event_position)
ggplot(df, aes(x = proportion_tree2, y = value)) +
  geom_point() +
  facet_grid(variable~tree_age, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# non reciprocal balanced
df <- subset(plot2_df, tree1_tree_shape == 'balanced')
df <- subset(df, tree2_tree_shape == 'balanced')
df <- subset(df, tree2_event_type == 'nonreciprocal')
df <- subset(df, tree2_event_position == 'close')
df$type = paste(df$tree2_event_type, df$tree2_event_position)
ggplot(df, aes(x = proportion_tree2, y = value)) +
  geom_point() +
  facet_grid(variable~tree_age, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Plot 3
plot3_df["event_asfactor"] <- as.factor(plot3_df$number_of_events)
ggplot(data = plot3_df, aes(x = event_asfactor, y = value)) + geom_boxplot(data = plot3_df, aes(fill = as.factor(tree_age))) + facet_wrap(~variable,scales = "free")
# attempt 2, reciprocal
df <- subset(plot3_df, tree1_tree_shape == 'balanced')
df <- subset(df, tree2_tree_shape == 'balanced')
df <- subset(df, tree2_event_type == 'reciprocal')
df <- subset(df, tree2_event_position == 'close')
df$type = paste(df$tree2_event_type, df$tree2_event_position)
ggplot(df, aes(x = event_asfactor, y = value)) +
  geom_boxplot() +
  facet_grid(variable~tree_age, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# attempt two, nonreciprocal
df <- subset(plot3_df, tree1_tree_shape == 'balanced')
df <- subset(df, tree2_tree_shape == 'balanced')
df <- subset(df, tree2_event_type == 'nonreciprocal')
df <- subset(df, tree2_event_position == 'close')
df$type = paste(df$tree2_event_type, df$tree2_event_position)
ggplot(df, aes(x = event_asfactor, y = value)) +
  geom_boxplot() +
  facet_grid(variable~tree_age, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Plot 4
ggplot(data=bs_df, aes(x="" ,y = value, fill = as.factor(tree_age))) + geom_boxplot(data = bs_df, aes(fill = as.factor(tree_age))) + facet_grid(as.factor(proportion_tree2)~variable)
ggplot(data=bs_df, aes(x=value)) + geom_histogram() + facet_grid(as.factor(proportion_tree2)~variable)
# attempt 2, reciprocal
df <- subset(bs_df, tree1_tree_shape == 'balanced')
df <- subset(df, tree2_tree_shape == 'balanced')
df <- subset(df, tree2_event_type == 'reciprocal')
df <- subset(df, tree2_event_position == 'close')
df$type = paste(df$tree2_event_type, df$tree2_event_position)
ggplot(df, aes(x = as.factor(proportion_tree2), y = value)) +
  geom_boxplot() +
  facet_grid(variable~tree_age, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  
ggplot(df, aes(x = value)) +
  geom_histogram() +
  facet_grid(variable~as.factor(proportion_tree2), scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  
# attempt 2, nonreciprocal
df <- subset(bs_df, tree1_tree_shape == 'balanced')
df <- subset(df, tree2_tree_shape == 'balanced')
df <- subset(df, tree2_event_type == 'nonreciprocal')
df <- subset(df, tree2_event_position == 'close')
df$type = paste(df$tree2_event_type, df$tree2_event_position)
ggplot(df, aes(x = as.factor(proportion_tree2), y = value)) +
  geom_boxplot() +
  facet_grid(variable~tree_age, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
# histogram?!?!?!?!!
ggplot(df, aes(x = value)) +
  geom_histogram() +
  facet_grid(variable~as.factor(proportion_tree2), scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))    

# no event plots
df <- subset(bs_df, tree1_tree_shape == 'balanced')
df <- subset(df, tree2_tree_shape == 'balanced')
df <- subset(df, tree2_event_type == 'none')
df <- subset(df, tree2_event_position == 'none')
df$type = paste(df$tree2_event_type, df$tree2_event_position)
ggplot(df, aes(x = as.factor(proportion_tree2), y = value)) +
  geom_boxplot() +
  facet_grid(variable~tree_age, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
# histogram?!?!?!?!!
ggplot(df, aes(x = value)) +
  geom_histogram() +
  facet_grid(variable~as.factor(proportion_tree2), scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Plots for thesis/seminar
# Plot 1: How do different events impact detection of recombination?
e = subset(plot1_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1 )
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree-splits proportion","pdm_difference" = "Summed element-wise difference","PHI_observed" = "PHI","prop_resolved_quartets" = "Proportion of resolved quartets", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets", "splittable_percentage" = "Modified splittable percentage")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
ggplot(e, aes(x = type, y = value)) +
  geom_boxplot() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_discrete(name = "Type of introgression event",
                   labels=c("none none" = "None", "reciprocal close" = "Reciprocal, \n Close", "reciprocal divergent" = "Reciprocal, \n Divergent", "reciprocal ancient" = "Reciprocal, \n Ancient",
                            "nonreciprocal close" = "Nonreciprocal, \n Close", "nonreciprocal divergent" = "Nonreciprocal, \n Divergent", "nonreciprocal ancient" = "Nonreciprocal, \n Ancient"),
                   limits=c("none none","reciprocal close","nonreciprocal close","reciprocal divergent","nonreciprocal divergent","reciprocal ancient","nonreciprocal ancient")) +
  ylab("Test statistic value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40))
# To get proper test size etc, save with the following dimensions: 4090 x 1938

# Plot 2: How does increasing the proportion of the recombinant sequence affect detection of treelikeness?
e = subset(plot2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree-splits proportion","pdm_difference" = "Summed element-wise difference","PHI_observed" = "PHI","prop_resolved_quartets" = "Proportion of resolved quartets", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets", "splittable_percentage" = "Modified splittable percentage")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
ggplot(e, aes(x = proportion_tree2, y = value)) +
  geom_point() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "Proportion of tree 2") +
  ylab("Test statistic value") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40))

# Plot 3: How does tree age affect detection of treelikeness?
e = subset(plot2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e$age = factor(e$tree_age)
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree-splits proportion","pdm_difference" = "Summed element-wise difference","PHI_observed" = "PHI","prop_resolved_quartets" = "Proportion of resolved quartets", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets", "splittable_percentage" = "Modified splittable percentage")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
ggplot(e, aes(x = proportion_tree2, y = value, color = age )) +
  geom_smooth() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "Proportion of tree 2") +
  ylab("Test statistic value") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40)) +
  scale_color_viridis(discrete = TRUE, option = "viridis")

# Plot 4: How does the number of events impact detection of tree likeness?
e = subset(plot3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree-splits proportion","pdm_difference" = "Summed element-wise difference","PHI_observed" = "PHI","prop_resolved_quartets" = "Proportion of resolved quartets", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets", "splittable_percentage" = "Modified splittable percentage")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
ggplot(e, aes(x = event_asfactor, y = value)) +
  geom_boxplot() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_discrete(name = "Number of introgression events") +
  ylab("Test statistic value") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40))
# To get proper test size etc, save with the following dimensions: 4090 x 1938

# Plot 5: How does reciprocity of events influence detection of treelikeness?
e = subset(plot3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e$group = factor(e$variable,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","splittable_percentage","pdm_difference","neighbour_net"))
facet_names <- list("neighbour_net" = "Tree-splits proportion","pdm_difference" = "Summed element-wise difference","PHI_observed" = "PHI","prop_resolved_quartets" = "Proportion of resolved quartets", 
                    "proportion_recombinant_triplets" = "Proportion of recombinant triplets", "splittable_percentage" = "Modified splittable percentage")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
ggplot(e, aes(x = number_of_events, y = value, colour = tree2_event_type)) +
  geom_smooth() +
  facet_wrap(~group,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_discrete(name = "Number of introgression events") +
  ylab("Test statistic value") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40)) + 
  scale_color_viridis(discrete = TRUE, option = "plasma")

# Plot 6: Are the results statistically significant?
e = subset(bs_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e$group = factor(e$variable,levels = c("PHI_p_value","X3Seq_p_value","likelihood_mapping_p_value","splittable_percentage_p_value","pdm_difference_p_value","neighbour_net_p_value"))
facet_names <- list("neighbour_net_p_value" = "Tree-splits \n proportion","pdm_difference_p_value" = "Summed element-wise \n difference","PHI_p_value" = "PHI","likelihood_mapping_p_value" = "Proportion of \n resolved quartets", 
                    "X3Seq_p_value" = "Proportion of \n recombinant triplets", "splittable_percentage_p_value" = "Modified splittable \n percentage")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~proportion_tree2,scales = "free_y", labeller = labeller(group = facet_labeller)) +
  scale_x_continuous(name = "P value") +
  ylab("Count") +
  theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 40), axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 30), strip.text = element_text(size = 40))

