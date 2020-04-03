# R code to import dataframes of test statistic and statistical test results and create plots to summarise the results
# Sourcing this file will open five dataframes and output a number of plots
# Final result is a number of plots displaying test statistic values under perturbation of various simulation factors



##### Step 1: Open packages #####
library(ggplot2)
library(gridExtra)
library(viridis)



##### Step 2: Specify file paths #####
# results_folder <- the folder where the result csvs will be placed (I use same results_folder in Parts 1 - 4.)
# plots_folder <- the folder where the plots will be stored
# maindir <- "treelikeness" repository location
# run_id <- if "run.id  = FALSE", program extracts run_id from input parameter file names 
#        <- otherwise, run_id will be set to whatever the user inputs here (e.g. "run_id = 'replicateAnalysis' ")

results_folder <- "/Users/caitlincherryh/Documents/Honours/Results/simulations_20200304/output/"
plots_folder <- "/Users/caitlincherryh/Documents/Honours/Results/simulations_20200304/plots/"
maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/"
run_id = FALSE



##### Step 3: Source function files #####
source(paste0(maindir,"code/func_process_data.R"))
# Extract run.id from the results folder name so the whole analysis has the same run.id
if (run_id == "FALSE"){
  run_id <- extract.run.id(results_folder)
}



##### Step 4: Open dataframes #####
# List all the files in the results_folder
melt_files <- list.files(results_folder)[grep("melted",list.files(results_folder))]
# Using the list of files in the results_folder, extract the melted csvs for each experiment (they are ready to plot!)
# p1:p3 are the dataframes for the test statistics calculated and estimated in experiments 1:3 respectively
ts1_df <- read.csv(paste0(results_folder,melt_files[grep("exp1_testStatistics",melt_files)]), stringsAsFactors = FALSE)
ts2_df <- read.csv(paste0(results_folder,melt_files[grep("exp2_testStatistics",melt_files)]), stringsAsFactors = FALSE)
ts3_df <- read.csv(paste0(results_folder,melt_files[grep("exp3_testStatistics",melt_files)]), stringsAsFactors = FALSE)
# bootstrap dataframe contains information about the p values (obtained for tree proportion using a parametric bootstrap)
# bs2:bs3 are the dataframes for the p values calculated and estimated in experiments 2:3 respectively
bs2_df <- read.csv(paste0(results_folder,melt_files[grep("exp2_p_value",melt_files)]), stringsAsFactors = FALSE)
bs3_df <- read.csv(paste0(results_folder,melt_files[grep("exp3_p_value",melt_files)]), stringsAsFactors = FALSE)



##### Step 5: Plot  #####
# Plot 1: How do different events impact detection of recombination?
e <- ts1_df # take the first experiment results and examine how different events impact detection of recombination
# Take only simulations that had a balanced tree 1 and tree 2
e = subset(e, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
# Fix the substitutions per site rate
e = subset(e, tree_age == 1 )
e = e[e$variable %in% c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"),]
# Combine event type and position into one column that can be used as an x-axis for the box plots
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e = e[e$type %in% c("none none","nonreciprocal close","reciprocal close"),]
# Create group to facet by
e$group <- factor(e$variable, levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q",
                                         "mode_delta_q","neighbour_net_trimmed"))
# Reorder variables so the facet grid comes out desired order - use new factor column
facet_names <- list("PHI_observed" = "PHI \n (PhiPack)","proportion_recombinant_triplets" = "Proportion of \n recombinant triplets \n (3SEQ)",
                    "prop_resolved_quartets" = "Proportion of resolved \n quartets \n (IQ-Tree)","neighbour_net_trimmed" = "Tree proportion \n (this paper)",
                    "mean_delta_q" = "Mean \u03B4q \n (\u03B4 plots)","mode_delta_q" = "Mode \u03B4q \n (\u03B4 plots)")
# Create a function to label facets
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
# Plot the events for each test statistic with both a fixed and a free y axis
# Fixed y plot:
p <- ggplot(e, aes(x = type, y = value)) +
  geom_boxplot(outlier.size = 3) +
  facet_wrap(~group, labeller = labeller(group = facet_labeller), ncol=3) +
  scale_x_discrete(name = "\n Type of introgression event \n",
                   labels=c("none none" = "None", "reciprocal close" = "Reciprocal, Close", "nonreciprocal close" = "Nonreciprocal, Close"),
                   limits=c("none none","reciprocal close","nonreciprocal close")) +
  ylab("\n Test statistic value \n") +
  theme(axis.title.x = element_text(size = 35), axis.title.y = element_text(size = 35),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 30), axis.text.y = element_text(size = 30), 
        strip.text = element_text(size = 30), strip.text.x = element_text(margin = margin(1,1,1,1, "cm")))
ggsave(filename = paste0(plots_folder,"exp1_differentEventTypes_fixedy.png"), plot = p, units = "in", width = 16, height = 20)

# Free y plot:
p <- ggplot(e, aes(x = type, y = value)) +
  geom_boxplot(outlier.size = 3) +
  facet_wrap(~group, scale = "free_y", labeller = labeller(group = facet_labeller), ncol=3) +
  scale_x_discrete(name = "\n Type of introgression event \n",
                   labels=c("none none" = "None", "reciprocal close" = "Reciprocal, Close", "nonreciprocal close" = "Nonreciprocal, Close"),
                   limits=c("none none","reciprocal close","nonreciprocal close")) +
  ylab("\n Test statistic value \n") +
  theme(axis.title.x = element_text(size = 35), axis.title.y = element_text(size = 35),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 30), axis.text.y = element_text(size = 30), 
        strip.text = element_text(size = 30), strip.text.x = element_text(margin = margin(1,1,1,1, "cm")))
ggsave(filename = paste0(plots_folder,"exp1_differentEventTypes_freey.png"), plot = p, units = "in", width = 17, height = 20)


############UPDATE THIS PLOT CAITLIN######################
# Plot 2: How does increasing the proportion of the recombinant sequence affect detection of treelikeness?
e = subset(ts3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e = e[e$variable %in% c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","pdm_difference","proportion_recombinant_triplets","neighbour_net_untrimmed",
                                       "neighbour_net_trimmed","prop_resolved_quartets","mean_delta_q","mode_delta_q"))
facet_names <- list("PHI_observed" = "PHI \n (PhiPack)","proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)",
                    "prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", "splittable_percentage" = "Distance ratio \n (this paper)",
                    "pdm_difference" = "Distance difference \n (this paper)","neighbour_net_untrimmed" = "Tree proportion \n (Untrimmed) \n (this paper)",
                    "neighbour_net_trimmed" = "Tree proportion \n (Trimmed) \n (this paper)", "mean_delta_q" = "Mean delta_q \n (delta plots)",
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

############UPDATE THIS PLOT CAITLIN######################
# Plot 3: How does tree age affect detection of treelikeness?
print("Plot 3")
e = subset(ts3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$age = factor(e$tree_age)
e = subset(e, variable == "PHI_observed" | variable == "splittable_percentage"  | variable == "pdm_difference" | variable == "proportion_recombinant_triplets" |
             variable == "neighbour_net_untrimmed" | variable == "neighbour_net_trimmed" | variable == "prop_resolved_quartets" | variable == "mean_delta_q" |
             variable == "mode_delta_q")
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","pdm_difference","proportion_recombinant_triplets","neighbour_net_untrimmed",
                                       "neighbour_net_trimmed","prop_resolved_quartets","mean_delta_q","mode_delta_q"))
facet_names <- list("PHI_observed" = "PHI \n (PhiPack)","proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)",
                    "prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", "splittable_percentage" = "Distance ratio \n (this paper)",
                    "pdm_difference" = "Distance difference \n (this paper)","neighbour_net_untrimmed" = "Tree proportion \n (Untrimmed) \n (this paper)",
                    "neighbour_net_trimmed" = "Tree proportion \n (Trimmed) \n (this paper)", "mean_delta_q" = "Mean delta_q \n (delta plots)",
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

############UPDATE THIS PLOT CAITLIN######################
# Plot 4: How does the number of events impact detection of tree likeness?
print("Plot 4")
e = subset(ts2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "reciprocal")
e$event_asfactor <- as.factor(e$number_of_events)
e = subset(e, variable == "PHI_observed" | variable == "splittable_percentage"  | variable == "pdm_difference" | variable == "proportion_recombinant_triplets" |
             variable == "neighbour_net_untrimmed" | variable == "neighbour_net_trimmed" | variable == "prop_resolved_quartets" | variable == "mean_delta_q" |
             variable == "mode_delta_q")
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","pdm_difference","proportion_recombinant_triplets","neighbour_net_untrimmed",
                                       "neighbour_net_trimmed","prop_resolved_quartets","mean_delta_q","mode_delta_q"))
facet_names <- list("PHI_observed" = "PHI \n (PhiPack)","proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)",
                    "prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", "splittable_percentage" = "Distance ratio \n (this paper)",
                    "pdm_difference" = "Distance difference \n (this paper)","neighbour_net_untrimmed" = "Tree proportion \n (Untrimmed) \n (this paper)",
                    "neighbour_net_trimmed" = "Tree proportion \n (Trimmed) \n (this paper)", "mean_delta_q" = "Mean delta_q \n (delta plots)",
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

############UPDATE THIS PLOT CAITLIN######################
# Plot 5: How does reciprocity of events influence detection of treelikeness?
print("Plot 5")
e = subset(ts3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, variable == "PHI_observed" | variable == "splittable_percentage"  | variable == "pdm_difference" | variable == "proportion_recombinant_triplets" |
             variable == "neighbour_net_untrimmed" | variable == "neighbour_net_trimmed" | variable == "prop_resolved_quartets" | variable == "mean_delta_q" |
             variable == "mode_delta_q")
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group = factor(e$variable,levels = c("PHI_observed","splittable_percentage","pdm_difference","proportion_recombinant_triplets","neighbour_net_untrimmed",
                                       "neighbour_net_trimmed","prop_resolved_quartets","mean_delta_q","mode_delta_q"))
facet_names <- list("PHI_observed" = "PHI \n (PhiPack)","proportion_recombinant_triplets" = "Proportion of recombinant triplets \n (3SEQ)",
                    "prop_resolved_quartets" = "Proportion of resolved quartets \n (IQ-Tree)", "splittable_percentage" = "Distance ratio \n (this paper)",
                    "pdm_difference" = "Distance difference \n (this paper)","neighbour_net_untrimmed" = "Tree proportion \n (Untrimmed) \n (this paper)",
                    "neighbour_net_trimmed" = "Tree proportion \n (Trimmed) \n (this paper)", "mean_delta_q" = "Mean delta_q \n (delta plots)",
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

############UPDATE THIS PLOT CAITLIN######################
# Plot 6: increasing proportion of introgressed DNA - are the results statistically significant?
# Histograms for p-values, for 0% - 50% introgression in 10% increments, for all test statistics
e = subset(bs3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e = subset(e, variable != "PHI_observed_p_value")
e = subset(e, variable != "num_recombinant_sequences_p_value")

e = subset(e, variable == "PHI_p_value" | variable == "splittable_percentage_p_value"  | variable == "pdm_difference_p_value" | variable == "X3Seq_p_value" |
             variable == "neighbour_net_trimmed_p_value" | variable == "neighbour_net_untrimmed_p_value" | variable == "likelihood_mapping_p_value" | variable == "mean_delta_q_p_value" |
             variable == "mode_delta_q_p_value")
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
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

############UPDATE THIS PLOT CAITLIN######################
# Plot 7: increasing proportion of introgressed DNA - are the results statistically significant?
# Line graph for p-values, for 0% - 50% introgression in 10% increments, for all test statistics
e = subset(bs3_df, tree1_tree_shape == 'balanced')
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


# Plot 8: increasing number of introgression events - are the results statistically significant?
# Histograms for p-values, for 0% - 50% introgression in 10% increments, for all test statistics
e = subset(bs2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "reciprocal")
e = e[e$variable %in% c("PHI_p_value","X3Seq_p_value","likelihood_mapping_p_value","mean_delta_q_p_value","mode_delta_q_p_value","neighbour_net_trimmed_p_value"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group = factor(e$variable,levels = c("PHI_p_value","X3Seq_p_value","likelihood_mapping_p_value",
                                       "mean_delta_q_p_value","mode_delta_q_p_value","neighbour_net_trimmed_p_value"))
facet_names <- list("PHI_p_value" = "PHI \n (PhiPack)","X3Seq_p_value" = "3SEQ",
                    "neighbour_net_trimmed_p_value" = "Tree \n proportion \n (this paper)","likelihood_mapping_p_value" = "Proportion of \n resolved quartets \n (IQ-Tree)",
                    "mean_delta_q_p_value" = "Mean \u03B4q \n (\u03B4 plots)","mode_delta_q_p_value" = "Mode \u03B4q \n (\u03B4 plots)")

facet_labeller <- function(variable){
  variable <- facet_names[variable]
}
p <- ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~number_of_events, labeller = labeller(group = facet_labeller)) +
  xlab("\n p  value \n") +
  ylab("\n Count \n") +
  theme(axis.text.x = element_text(size = 55), axis.title.x = element_text(size = 90), axis.title.y = element_text(size = 90),
        axis.text.y = element_text(size = 55), strip.text = element_text(size = 60), strip.text.x = element_text(margin = margin(1,1,1,1, "cm")),
        strip.text.y = element_text(margin = margin(1,1,1,1, "cm"))) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1)) + 
  scale_y_continuous(breaks = seq(0,10,2), labels = seq(0,10,2))
ggsave(filename = paste0(plots_folder,"exp2_IncNumEvents_StatisticalSignificance_hists.png"), plot = p, units = "in", width = 70, height =53, limitsize = FALSE)


# Plot 9: increasing number of introgression events - are the results statistically significant?
# Line graph for p-values, for 0% - 50% introgression in 10% increments, for all test statistics
e = subset(bs2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == 1)
e = subset(e, tree2_event_type != "reciprocal")
# extract values from the dataframe for the new plot
PHI_p_value <- c(nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 0,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 1,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 2,]),
                 nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 3,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 4,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 5,]),
                 nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 6,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 7,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 8,]))
X3SEQ_p_value <- c(nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 0,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 1,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 2,]),
                   nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 3,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 4,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 5,]),
                   nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 6,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 7,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 8,]))
likelihood_mapping_p_value <- c(nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 0,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 1,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 2,]),
                                nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 3,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 4,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 5,]),
                                nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 6,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 7,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 8,]))
neighbour_net_trimmed_p_value <- c(nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 0,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 1,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 2,]),
                                   nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 3,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 4,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 5,]),
                                   nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 6,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 7,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 8,]))
mean_delta_q_p_value <- c(nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 0,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 1,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 2,]),
                          nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 3,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 4,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 5,]),
                          nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 6,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 7,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 8,]))
mode_delta_q_p_value <- c(nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 0,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 1,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 2,]),
                          nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 3,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 4,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 5,]),
                          nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 6,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 7,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 8,]))
# make a new dataframe for the new plot
value <- c(PHI_p_value,X3SEQ_p_value,likelihood_mapping_p_value,mean_delta_q_p_value,mode_delta_q_p_value,neighbour_net_trimmed_p_value)
# get value as a percentage by dividing by 10 (number of sequences per event)
value <- value/10*100
ts <- c(rep("PHI_p_value",9), rep("X3SEQ_p_value",9), rep("likelihood_mapping_p_value",9), rep("mean_delta_q_p_value",9), rep("mode_delta_q_p_value",9),rep("neighbour_net_trimmed_p_value",9))
num_events <- c(rep(seq(0,8,1),6))
f <- data.frame(num_events,ts,value,stringsAsFactors = FALSE)
f$variable <- factor(ts, levels = c("PHI_p_value", "X3SEQ_p_value", "likelihood_mapping_p_value", "mean_delta_q_p_value", "mode_delta_q_p_value",
                                    "neighbour_net_trimmed_p_value"))
facet_names <- list("PHI_p_value" = "PHI \n (PhiPack)", "likelihood_mapping_p_value" = "Proportion of \n resolved quartets \n (IQ-Tree)",
                    "neighbour_net_trimmed_p_value" = "Tree \n proportion \n (this paper)", "X3SEQ_p_value" = "3SEQ",
                    "mean_delta_q_p_value" = "Mean  \u03B4q \n (\u03B4 plots)","mode_delta_q_p_value" = "Mode  \u03B4q \n (\u03B4 plots)")
facet_labeller <- function(variable){
  variable <- facet_names[variable]
}

p <- ggplot(f, aes(x = num_events, y = value)) +
  geom_line(size=3) +
  facet_wrap(~variable, labeller = labeller(variable = facet_labeller), nrow = 2, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n", labels = seq(0,8,1), breaks = seq(0,8,1)) +
  ylab("\n Percent of simulations that reject the null hypothesis \n (p-value < 0.05) \n") +
  theme(axis.text.x = element_text(size = 40), axis.title.x = element_text(size = 50), axis.title.y = element_text(size = 50),
        axis.text.y = element_text(size = 40), strip.text = element_text(size = 50), strip.text.x = element_text(margin = margin(1,0,0.5,0, "cm")),
        strip.text.y = element_text(margin = margin(1,0,0.5,0, "cm")), legend.text = element_text(size = 40),
        legend.title = element_text(size = 50), legend.key.width = unit(5,"cm"), legend.key.height = unit(5, "cm"),
        panel.background = element_rect(fill="white"),panel.grid.major = element_line(colour = "#999999"),panel.grid.minor = element_line(colour = "grey78")) +
  scale_y_continuous(labels = seq(0,100,10), breaks = seq(0,100,10), minor_breaks = seq(0,100,5), limits = c(0,100)) + 
  geom_hline(aes(yintercept = 5, colour = "red"), linetype = "dashed", size = 1.5) + 
  scale_colour_manual("Ideal false\npositive rate\n", values="red", labels = "5% when\n\u03b1 = 0.05")
ggsave(filename = paste0(plots_folder,"exp2_IncNumEvents_StatisticalSignificance_line.png"), plot = p, units = "in", width = 40, height = 46.8)
