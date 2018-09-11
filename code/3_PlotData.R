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

# Open dataframes
plot1_df <- read.csv(paste0(input_folder,"plot1_testStatistics_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)
plot2_df <- read.csv(paste0(input_folder,"plot2_testStatistics_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)
plot3_df <- read.csv(paste0(input_folder,"plot3_testStatistics_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)
plot4_df <- read.csv(paste0(input_folder,"plot4_testStatistics_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)
bs_df <-  read.csv(paste0(input_folder,"plot4_p_value_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)


# Plot 1:
ggplot(data=plot1_df,aes(x=tree2_event_position,y=value)) + geom_boxplot(data=plot1_df,aes(fill = tree2_event_type)) + facet_wrap(~variable) + ylim(0,2)

# Plot 2
ggplot(data=plot2_df,aes(x=proportion_tree2,y=value)) + geom_point(data=plot2_df, aes(colour = as.factor(tree_age))) + facet_wrap(~variable) + geom_smooth(method = "lm", data = plot2_df, aes(x=proportion_tree2,y=value, colour = as.factor(tree_age)))
ggplot(data=plot2_df,aes(x=proportion_tree2,y=value)) + facet_wrap(~variable) + geom_smooth(method = "lm", data = plot2_df, aes(x=proportion_tree2,y=value, colour = as.factor(tree_age)))

# Plot 3
plot3_df["event_asfactor"] <- as.factor(plot3_df$number_of_events)
ggplot(data = plot3_df, aes(x = event_asfactor, y = value)) + geom_boxplot(data = plot3_df, aes(fill = as.factor(tree_age))) + facet_wrap(~variable,scales = "free")

# Plot 4
ggplot(data=bs_df, aes(x="" ,y = value, fill = as.factor(tree_age))) + geom_boxplot(data = bs_df, aes(fill = as.factor(tree_age))) + facet_grid(as.factor(proportion_tree2)~variable)

       