# R code to import dataframes of test statistic and statistical test results and create plots to summarise the results
# Sourcing this file will open five dataframes and output a number of plots
# Final result is a number of plots displaying test statistic values under perturbation of various simulation factors

# completed plots: 1,2

##### Step 1: Open packages #####
library(ggplot2)
library(ggpmisc)
library(patchwork)




##### Step 2: Uncomment and set the file paths for output folders, executables, and the identifying name for this run #####
# results_folder <- the folder where the result csvs will be placed (I use same results_folder in Parts 1 - 4.)
# plots_folder <- the folder where the plots will be stored
# maindir <- "treelikeness" repository location
# run_id <- if "run.id  = FALSE", program extracts run_id from input parameter file names 
#        <- otherwise, run_id will be set to whatever the user inputs here (e.g. "run_id = 'replicateAnalysis' ")

# results_folder <- ""
# plots_folder <- ""
# maindir <- ""
# run_id = FALSE
# tree_length = 0.5  # choose a value for total tree depth: 0.05, 0.10, 0.5, or 1
                     # this value will be used as tree depth in all plots EXCEPT the plots that contain all 4 tree depths for comparison


#__________________________________________Caitlin's paths (delete these if you're not Caitlin)______________________________________
results_folder <- "/Users/caitlincherryh/Documents/Honours/Results/simulations_20200304/output/"
plots_folder <- "/Users/caitlincherryh/Documents/Honours/Results/simulations_20200304/plots/"
maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/"
run_id = FALSE
tree_length = 0.5
#____________________________________________________________________________________________________________________________________



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
# Old facet_names and labeller in case it's needed later on:
# facet_names <- list("PHI_observed" = "PHI \n (PhiPack)","proportion_recombinant_triplets" = "Proportion of \n recombinant triplets \n (3SEQ)",
#                     "prop_resolved_quartets" = "Proportion of resolved \n quartets \n (IQ-Tree)","neighbour_net_trimmed" = "Tree proportion \n (this paper)",
#                     "mean_delta_q" = "Mean \u03B4q \n (\u03B4 plots)","mode_delta_q" = "Mode \u03B4q \n (\u03B4 plots)")
# facet_labeller <- function(variable){
#   variable <- facet_names[variable]
# }
# Panel colour elements: panel.background = element_rect(fill="white"), panel.grid.major = element_line(colour = "#999999"),panel.grid.minor = element_line(colour = "grey78")


# Plot 1: How do different events impact detection of recombination?
print("Plot 1")
e <- ts1_df # take the first experiment results and examine how different events impact detection of recombination
# Take only simulations that had a balanced tree 1 and tree 2
e = subset(e, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
# Fix the substitutions per site rate
e = subset(e, tree_age == tree_length)
e = e[e$variable %in% c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"),]
# Combine event type and position into one column that can be used as an x-axis for the box plots
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e = e[e$type %in% c("none none","nonreciprocal close","reciprocal close"),]
# Create group to facet by
e$group <- factor(e$variable, levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Proportion of recombinant triplets","(3SEQ)")),
                             expression(atop("Proportion of resolved quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

# Plot the events for each test statistic with both a fixed and a free y axis
# Fixed y plot:
p <- ggplot(e, aes(x = type, y = value)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~group, labeller = label_parsed, ncol=3) +
  scale_x_discrete(name = "\n Type of introgression event \n",
                   labels=c("none none" = "None", "reciprocal close" = "Reciprocal, Close", "nonreciprocal close" = "Nonreciprocal, Close"),
                   limits=c("none none","reciprocal close","nonreciprocal close")) +
  ylab("\n Test statistic value \n") + theme_bw() + 
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 10), 
        strip.text = element_text(size = 9), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp1_differentEventTypes_fixedy.png"), plot = p, units = "in", width = 7, height = 7)


cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp1_differentEventTypes_fixedy.pdf"), height = 7, width = 7, fallback_resolution = 300)
ggplot(e, aes(x = type, y = value)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~group, labeller = label_parsed, ncol=3) +
  scale_x_discrete(name = "\n Type of introgression event \n",
                   labels=c("none none" = "None", "reciprocal close" = "Reciprocal, Close", "nonreciprocal close" = "Nonreciprocal, Close"),
                   limits=c("none none","reciprocal close","nonreciprocal close")) +
  ylab("\n Test statistic value \n") + theme_bw() + 
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 10), 
        strip.text = element_text(size = 9), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")))
dev.off()

# Free y plot:
p <- ggplot(e, aes(x = type, y = value)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~group, scale = "free_y", labeller = label_parsed, ncol=3) +
  scale_x_discrete(name = "\n Type of introgression event \n",
                   labels=c("none none" = "None", "reciprocal close" = "Reciprocal, Close", "nonreciprocal close" = "Nonreciprocal, Close"),
                   limits=c("none none","reciprocal close","nonreciprocal close")) +
  ylab("\n Test statistic value \n") + theme_bw() + 
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 10), 
        strip.text = element_text(size = 9), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp1_differentEventTypes_freey.png"), plot = p, units = "in", height = 8, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp1_differentEventTypes_freey.pdf"), height = 8, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = type, y = value)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~group, scale = "free_y", labeller = label_parsed, ncol=3) +
  scale_x_discrete(name = "\n Type of introgression event \n",
                   labels=c("none none" = "None", "reciprocal close" = "Reciprocal, Close", "nonreciprocal close" = "Nonreciprocal, Close"),
                   limits=c("none none","reciprocal close","nonreciprocal close")) +
  ylab("\n Test statistic value \n") + theme_bw() + 
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 10), 
        strip.text = element_text(size = 9), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")))
dev.off()




# Plot 2: How does increasing the proportion of the recombinant sequence affect detection of treelikeness?
print("Plot 2")
e = subset(ts3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == tree_length)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e = e[e$variable %in% c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group <- factor(e$variable, levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Proportion of recombinant triplets","(3SEQ)")),
                             expression(atop("Proportion of resolved quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

# create a vector of the r^2 values from each variable
r2_raw <- c()
var_list <- c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed")
for (var in var_list){
  # get only entries with this variable
  var_df <- e[(e$variable == var),]
  # insert into regression model
  m <- lm(value ~ proportion_tree2, var_df)
  r2_var <- summary(m)$r.squared
  r2_raw <- c(r2_raw, r2_var)
}

# Quick function to reformat r^2 values the way I like
reformat.r2 <- function(number){
  number <- round(number, digits = 3)
}

# Apply function to r^2 values
r2_plot <- unlist(lapply(r2_raw,reformat.r2))
replace_inds <- which(r2_plot == 0)
r2_plot <- as.character(r2_plot)
r2_plot[replace_inds] <- "0.000" # replace the 0 with 0.000 (as R^2 was calculated to 3dp)
r2_plot <- paste0(" = ", r2_plot)

# Make a dataframe of variables r^2 values, vectors for plot placement and equation for the plot
e_eq <- data.frame(variable = var_list, x_pos = rep(0.0, 6) , y1_pos = rep(1.1,6), y2_pos = c(0.525, 0.0875, 1.01, 0.4, 0.5, 1),
                   r2_x_pos = rep(0.035, 6), rsq = r2_plot, rsquare = "R^2 ",
                   group = factor(var_list,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), 
                                  ordered = TRUE, 
                                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Proportion of recombinant triplets","(3SEQ)")),
                                             expression(atop("Proportion of resolved quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) ) ) 

# Plot points with regression line and r^2 value
p <- ggplot(e, aes(x = proportion_tree2, y = value)) +
  geom_point(colour = "gray55", alpha = 0.2) +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
  facet_wrap(~group, labeller = label_parsed, ncol=3) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) + 
  geom_text(data = e_eq, aes(x = x_pos, y = y1_pos, label = rsquare), parse = TRUE, hjust = 0, size = 3) + 
  geom_text(data = e_eq, aes(x = r2_x_pos, y = y1_pos, label = rsq), parse = FALSE, hjust = 0, size = 3)
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_increasingProportionTree2_regression_fixedy.png"), plot = p, units = "in", height = 6.6, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_increasingProportionTree2_regression_fixedy.pdf"), height = 6.57, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = proportion_tree2, y = value)) +
  geom_point(colour = "gray55", alpha = 0.2) +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
  facet_wrap(~group, labeller = label_parsed, ncol=3) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) + 
  geom_text(data = e_eq, aes(x = x_pos, y = y1_pos, label = rsquare), parse = TRUE, hjust = 0, size = 3) + 
  geom_text(data = e_eq, aes(x = r2_x_pos, y = y1_pos, label = rsq), parse = FALSE, hjust = 0, size = 3)
dev.off()

p <- ggplot(e, aes(x = proportion_tree2, y = value)) +
  geom_point(colour = "gray55", alpha = 0.2) +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
  facet_wrap(~group, scale = "free_y", labeller = label_parsed, ncol=3) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) + 
  geom_text(data = e_eq, aes(x = x_pos, y = y2_pos, label = rsquare), parse = TRUE, hjust = 0, size = 3) + 
  geom_text(data = e_eq, aes(x = r2_x_pos, y = y2_pos, label = rsq), parse = FALSE, hjust = 0, size = 3)
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_increasingProportionTree2_regression_freey.png"), plot = p, units = "in", height = 6.57, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_increasingProportionTree2_regression_freey.pdf"), height = 6.57, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = proportion_tree2, y = value)) +
  geom_point(colour = "gray55", alpha = 0.2) +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
  facet_wrap(~group, scale = "free_y", labeller = label_parsed, ncol=3) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) + 
  geom_text(data = e_eq, aes(x = x_pos, y = y2_pos, label = rsquare), parse = TRUE, hjust = 0, size = 3) + 
  geom_text(data = e_eq, aes(x = r2_x_pos, y = y2_pos, label = rsq), parse = FALSE, hjust = 0, size = 3)
dev.off()

# Continuation of plot 2, but includes all four tree lengths
e = subset(ts3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$age = factor(e$tree_age, ordered = TRUE)
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e = e[e$variable %in% c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group <- factor(e$variable, levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Prop. recomb. trip.","(3SEQ)")),
                             expression(atop("Prop. res. quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

r2_raw <- c()
var_list <- c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed")
eq_var_list <- rep(var_list, each = 4)
age_list <- c(0.05, 0.1, 0.5, 1)
for (var in var_list){
  var_df <- e[(e$variable == var),]
  for (t_age in age_list){
    # get only entries with this variable
    # insert into regression model
    age_df <- var_df[(var_df$age == t_age),]
    m <- lm(value ~ proportion_tree2, age_df)
    r2_var <- summary(m)$r.squared
    print(paste0("Variable: ",var," - tree depth: ",t_age," - R^2 : ",round(r2_var, digits = 3)))
    r2_raw <- c(r2_raw, r2_var) 
  }
}

# Quick function to reformat r^2 values the way I like
reformat.r2 <- function(number){
  number <- round(number, digits = 3)
}

# Apply function to r^2 values
r2_plot <- unlist(lapply(r2_raw,reformat.r2))
replace_inds <- which(r2_plot == 0)
r2_plot <- as.character(r2_plot)
r2_plot[replace_inds] <- "0.000" # replace the 0 with 0.000 (as R^2 was calculated to 3dp)
r2_plot <- paste0(" = ", r2_plot)

# Make a dataframe of variables r^2 values, vectors for plot placement and equation for the plot
e_eq <- data.frame(variable = eq_var_list, x_pos = rep(0.0, 24) , y1_pos = rep(1.1,24), y2_pos = rep(c(0.6, 0.11, 1.05, 0.45, 0.5, 1), each = 4), age = rep(c(0.05, 0.1, 0.5, 1.0), 6),
                   r2_x_pos = rep(0.035, 6), rsq = r2_plot, rsquare = "R^2 ",
                   group = factor(eq_var_list,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), 
                                  ordered = TRUE, 
                                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Prop. recomb. trip.","(3SEQ)")),
                                             expression(atop("Prop. res. quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) ) ) 
                                  
p <- ggplot(e, aes(x = proportion_tree2, y = value)) +
  geom_point(colour = "gray55", alpha = 0.2) +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
  facet_grid(group~age, scale = "free_y", labeller = label_parsed) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) + 
  geom_text(data = e_eq, aes(x = x_pos, y = y2_pos, label = rsquare), parse = TRUE, hjust = 0, size = 3) + 
  geom_text(data = e_eq, aes(x = r2_x_pos, y = y2_pos, label = rsq), parse = FALSE, hjust = 0, size = 3)
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_increasingProportionTree2_allTreeDepths_regression_freey.png"), plot = p, units = "in", height = 9, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_increasingProportionTree2_allTreeDepths_regression_freey.pdf"), height = 9, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = proportion_tree2, y = value)) +
  geom_point(colour = "gray55", alpha = 0.2) +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
  facet_grid(group~age, scale = "free_y", labeller = label_parsed) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) + 
  geom_text(data = e_eq, aes(x = x_pos, y = y2_pos, label = rsquare), parse = TRUE, hjust = 0, size = 3) + 
  geom_text(data = e_eq, aes(x = r2_x_pos, y = y2_pos, label = rsq), parse = FALSE, hjust = 0, size = 3)
dev.off()
  
patchwork <- p
cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_increasingProportionTree2_allTreeDepths_regression_freey_withTitle.pdf"), height = 9, width = 9, fallback_resolution = 300)
patchwork + 
  plot_annotation(title = "Total tree depth (substitutions per site)",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(family = "")))
dev.off()

# Continuation of plot 2, but includes both event types
e = subset(ts3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_event_type != "none")
e = subset(e, tree_age == tree_length)
e$tree2_event_type = factor(e$tree2_event_type, ordered = TRUE)
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e = e[e$variable %in% c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group <- factor(e$variable, levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Prop. recomb. trip.","(3SEQ)")),
                             expression(atop("Prop. res. quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

r2_raw <- c()
var_list <- c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed")
eq_var_list <- rep(var_list, each = 2)
event_list <- c("nonreciprocal","reciprocal")
for (var in var_list){
  var_df <- e[(e$variable == var),]
  for (ev in event_list){
    # get only entries with this variable
    # insert into regression model
    age_df <- var_df[(var_df$tree2_event_type == ev),]
    m <- lm(value ~ proportion_tree2, age_df)
    r2_var <- summary(m)$r.squared
    print(paste0("Variable: ",var," - event type: ",ev," - R^2 : ", round(r2_var, digits = 3) ))
    r2_raw <- c(r2_raw, r2_var) 
  }
}

# Quick function to reformat r^2 values the way I like
reformat.r2 <- function(number){
  number <- round(number, digits = 3)
}

# Apply function to r^2 values
r2_plot <- unlist(lapply(r2_raw,reformat.r2))
replace_inds <- which(r2_plot == 0)
r2_plot <- as.character(r2_plot)
r2_plot[replace_inds] <- "0.000" # replace the 0 with 0.000 (as R^2 was calculated to 3dp)
r2_plot <- paste0(" = ", r2_plot)


# Make a dataframe of variables r^2 values, vectors for plot placement and equation for the plot
e_eq <- data.frame(variable = eq_var_list, x_pos = rep(0.0, 12) , y1_pos = rep(1.1,12), 
                   y2_pos = rep(c(0.3875, 0.0875, 1.02, 0.28, 0.12, 1.05), each = 2),
                   r2_x_pos = rep(0.05, 6), rsq = r2_plot, rsquare = "R^2 ", tree2_event_type = rep(c("nonreciprocal","reciprocal"),6),
                   fac_event_type = factor(rep(event_list,6), ordered = TRUE, labels = c("Nonreciprocal","Reciprocal")),
                   group = factor(eq_var_list,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), 
                                  ordered = TRUE, 
                                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Prop. recomb. trip.","(3SEQ)")),
                                             expression(atop("Prop. res. quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) ) ) 
# Change e$tree2_event_type to a new factor with pretty formatting for output
e$fac_event_type <- factor(e$tree2_event_type, ordered = TRUE, labels = c("Nonreciprocal","Reciprocal"))

p <- ggplot(e, aes(x = proportion_tree2, y = value)) +
  geom_point(colour = "gray55", alpha = 0.2) +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
  facet_grid(group ~ fac_event_type, scale = "free_y", labeller = label_parsed) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  geom_text(data = e_eq, aes(x = r2_x_pos, y = y2_pos, label = rsq), parse = FALSE, hjust = 0, size = 3) +
  geom_text(data = e_eq, aes(x = x_pos, y = y2_pos, label = rsquare), parse = TRUE, hjust = 0, size = 3) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_increasingProportionTree2_RecipNonrecip_regression_freey.png"), plot = p, units = "in", height = 7, width = 6)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_increasingProportionTree2_RecipNonrecip_regression_freey.pdf"), height = 7, width = 6, fallback_resolution = 300)
ggplot(e, aes(x = proportion_tree2, y = value)) +
  geom_point(colour = "gray55", alpha = 0.2) +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
  facet_grid(group ~ fac_event_type, scale = "free_y", labeller = label_parsed) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  geom_text(data = e_eq, aes(x = r2_x_pos, y = y2_pos, label = rsq), parse = FALSE, hjust = 0, size = 3) +
  geom_text(data = e_eq, aes(x = x_pos, y = y2_pos, label = rsquare), parse = TRUE, hjust = 0, size = 3) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")))
dev.off()

patchwork <- p
cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_increasingProportionTree2_RecipNonrecip_regression_freey_withTitle.pdf"), height = 8, width = 6, fallback_resolution = 300)
patchwork + 
  plot_annotation(title = "Type of introgression event",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(family = "")))
dev.off()



# Plot 3: How does tree age affect detection of treelikeness?
print("Plot 3")
e = subset(ts3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e$age = factor(e$tree_age, ordered = TRUE)
e = e[e$variable %in% c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group <- factor(e$variable, levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Proportion of recombinant triplets","(3SEQ)")),
                             expression(atop("Proportion of resolved quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

p <- ggplot(e, aes(x = proportion_tree2, y = value, color = age )) +
  geom_smooth(size = 0.5, aes(linetype = age), method = "lm") +
  facet_wrap(~group,scales = "free_y", labeller = label_parsed, nrow = 3, ncol = 3) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8), legend.background = element_rect(linetype = 1, size = 0.2, colour = 1)) + 
  scale_color_manual(values = c("#a6611a","#dfc27d","#80cdc1","#018571"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00")) +
  scale_linetype_manual(values = c("longdash","longdash","solid","solid"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00"))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_treeAgeWithIncreasingTree2_lm_freey.png"),plot = p, units = "in", height = 6.6, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_treeAgeWithIncreasingTree2_lm_freey.pdf"), height = 6.6, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = proportion_tree2, y = value, color = age )) +
  geom_smooth(size = 0.5, aes(linetype = age), method = "lm") +
  facet_wrap(~group,scales = "free_y", labeller = label_parsed, nrow = 3, ncol = 3) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8), legend.background = element_rect(linetype = 1, size = 0.2, colour = 1)) + 
  scale_color_manual(values = c("#a6611a","#dfc27d","#80cdc1","#018571"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00")) +
  scale_linetype_manual(values = c("longdash","longdash","solid","solid"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00"))
dev.off()

p <- ggplot(e, aes(x = proportion_tree2, y = value, color = age )) +
  geom_smooth(size = 0.5, aes(linetype = age), method = "gam") +
  facet_wrap(~group,scales = "free_y", labeller = label_parsed, nrow = 3, ncol = 3) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8), legend.background = element_rect(linetype = 1, size = 0.2, colour = 1)) + 
  scale_color_manual(values = c("#a6611a","#dfc27d","#80cdc1","#018571"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00")) +
  scale_linetype_manual(values = c("longdash","longdash","solid","solid"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00"))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_treeAgeWithIncreasingTree2_gam_freey.png"),plot = p, units = "in", height = 6.6, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_treeAgeWithIncreasingTree2_gam_freey.pdf"), height = 6.6, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = proportion_tree2, y = value, color = age )) +
  geom_smooth(size = 0.5, aes(linetype = age), method = "gam") +
  facet_wrap(~group,scales = "free_y", labeller = label_parsed, nrow = 3, ncol = 3) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8), legend.background = element_rect(linetype = 1, size = 0.2, colour = 1)) + 
  scale_color_manual(values = c("#a6611a","#dfc27d","#80cdc1","#018571"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00")) +
  scale_linetype_manual(values = c("longdash","longdash","solid","solid"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00"))
dev.off()

p <- ggplot(e, aes(x = proportion_tree2, y = value, color = age )) +
  geom_smooth(size = 0.5, aes(linetype = age), method = "gam") +
  facet_wrap(~group, labeller = label_parsed, nrow = 3, ncol = 3) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8), legend.background = element_rect(linetype = 1, size = 0.2, colour = 1)) + 
  scale_color_manual(values = c("#a6611a","#dfc27d","#80cdc1","#018571"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00")) +
  scale_linetype_manual(values = c("longdash","longdash","solid","solid"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00"))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_treeAgeWithIncreasingTree2_gam_fixedy.png"),plot = p, units = "in", height = 6.6, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_treeAgeWithIncreasingTree2_gam_fixedy.pdf"), height = 6.6, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = proportion_tree2, y = value, color = age )) +
  geom_smooth(size = 0.5, aes(linetype = age), method = "gam") +
  facet_wrap(~group, labeller = label_parsed, nrow = 3, ncol = 3) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8), legend.background = element_rect(linetype = 1, size = 0.2, colour = 1)) + 
  scale_color_manual(values = c("#a6611a","#dfc27d","#80cdc1","#018571"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00")) +
  scale_linetype_manual(values = c("longdash","longdash","solid","solid"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00"))
dev.off()

e = subset(ts2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree2_event_type != "reciprocal")
e = subset(e, tree2_event_type != "none")
e$age = factor(e$tree_age, ordered = TRUE)
e = e[e$variable %in% c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group <- factor(e$variable, levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Proportion of recombinant triplets","(3SEQ)")),
                             expression(atop("Proportion of resolved quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

p <- ggplot(e, aes(x = number_of_events, y = value, color = age )) +
  geom_smooth(size = 0.5, aes(linetype = age), method = "loess") +
  facet_wrap(~group,scales = "free_y", labeller = label_parsed, nrow = 3, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n", breaks = c(0:8), labels = c(0:8)) +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 9), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8), legend.background = element_rect(linetype = 1, size = 0.2, colour = 1)) + 
  scale_color_manual(values = c("#a6611a","#dfc27d","#80cdc1","#018571"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00")) +
  scale_linetype_manual(values = c("longdash","longdash","solid","solid"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00"))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_treeAgeWithIncreasingNumEvents_loess_freey.png"),plot = p, units = "in", height = 6.6, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_treeAgeWithIncreasingNumEvents_loess_freey.pdf"), height = 6.6, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = number_of_events, y = value, color = age )) +
  geom_smooth(size = 0.5, aes(linetype = age), method = "loess") +
  facet_wrap(~group,scales = "free_y", labeller = label_parsed, nrow = 3, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n", breaks = c(0:8), labels = c(0:8)) +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 9), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8), legend.background = element_rect(linetype = 1, size = 0.2, colour = 1)) + 
  scale_color_manual(values = c("#a6611a","#dfc27d","#80cdc1","#018571"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00")) +
  scale_linetype_manual(values = c("longdash","longdash","solid","solid"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00"))
dev.off()

p <- ggplot(e, aes(x = number_of_events, y = value, color = age )) +
  geom_smooth(size = 0.5, aes(linetype = age), method = "loess") +
  facet_wrap(~group, labeller = label_parsed, nrow = 3, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n", breaks = c(0:8), labels = c(0:8)) +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8), legend.background = element_rect(linetype = 1, size = 0.2, colour = 1)) + 
  scale_color_manual(values = c("#a6611a","#dfc27d","#80cdc1","#018571"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00")) +
  scale_linetype_manual(values = c("longdash","longdash","solid","solid"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00"))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_treeAgeWithIncreasingNumEvents_loess_fixedy.png"),plot = p, units = "in", height = 6.6, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_treeAgeWithIncreasingNumEvents_loess_fixedy.pdf"), height = 6.6, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = number_of_events, y = value, color = age )) +
  geom_smooth(size = 0.5, aes(linetype = age), method = "loess") +
  facet_wrap(~group, labeller = label_parsed, nrow = 3, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n", breaks = c(0:8), labels = c(0:8)) +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        legend.text = element_text(size = 8), legend.title = element_text(size = 8), legend.background = element_rect(linetype = 1, size = 0.2, colour = 1)) + 
  scale_color_manual(values = c("#a6611a","#dfc27d","#80cdc1","#018571"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00")) +
  scale_linetype_manual(values = c("longdash","longdash","solid","solid"), name = expression(atop("Tree depth","(substitutions/site)")), labels = c("0.05","0.10","0.50","1.00"))
dev.off()




# Plot 4: How does the number of events impact detection of tree likeness?
print("Plot 4")
e = subset(ts2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == tree_length)
e = subset(e, tree2_event_type != "reciprocal")
e$event_asfactor <- as.factor(e$number_of_events)
e = e[e$variable %in% c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group <- factor(e$variable, levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Proportion of recombinant triplets","(3SEQ)")),
                             expression(atop("Proportion of resolved quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

p <- ggplot(e, aes(x = event_asfactor, y = value)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5, outlier.colour = "gray55", lwd = 0.4) +
  facet_wrap(~group, scales = "free_y", labeller = label_parsed, ncol = 3) +
  scale_x_discrete(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) 
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_numberOfEvents_freey.png"), plot = p, units = "in", height = 6.6, width = 9)


cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_numberOfEvents_freey.pdf"), height = 6.6, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = event_asfactor, y = value)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5, outlier.colour = "gray55", lwd = 0.4) +
  facet_wrap(~group, scales = "free_y", labeller = label_parsed, ncol = 3) +
  scale_x_discrete(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) 
dev.off()

p <- ggplot(e, aes(x = event_asfactor, y = value)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5, outlier.colour = "gray55", lwd = 0.4) +
  facet_wrap(~group, labeller = label_parsed, ncol = 3) +
  scale_x_discrete(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) 
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_numberOfEvents_fixedy.png"), plot = p, units = "in", height = 6.6, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_numberOfEvents_fixedy.pdf"), height = 6.6, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = event_asfactor, y = value)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5, outlier.colour = "gray55", lwd = 0.4) +
  facet_wrap(~group, labeller = label_parsed, ncol = 3) +
  scale_x_discrete(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) 
dev.off()

# Continuation of plot 4, but includes all four tree lengths
e = subset(ts2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_event_type != "reciprocal")
e$age = factor(e$tree_age, ordered = TRUE)
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e = e[e$variable %in% c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group <- factor(e$variable, levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Prop. recomb. trip.","(3SEQ)")),
                             expression(atop("Prop. res. quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

r2_raw <- c()
var_list <- c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed")
eq_var_list <- rep(var_list, each = 4)
age_list <- c(0.05, 0.1, 0.5, 1)
for (var in var_list){
  var_df <- e[(e$variable == var),]
  for (t_age in age_list){
    # get only entries with this variable
    # insert into regression model
    age_df <- var_df[(var_df$age == t_age),]
    age_df$age <- as.numeric(age_df$age)
    age_df$number_of_events <- as.numeric(age_df$number_of_events)
    m <- lm(value ~ number_of_events, age_df)
    r2_var <- summary(m)$r.squared
    print(paste0("Variable: ",var," - tree depth: ",t_age," - R^2 : ",round(r2_var, digits = 3)))
    r2_raw <- c(r2_raw, r2_var) 
  }
}

# Quick function to reformat r^2 values the way I like
reformat.r2 <- function(number){
  number <- round(number, digits = 3)
}

# Apply function to r^2 values
r2_plot <- unlist(lapply(r2_raw,reformat.r2))
replace_inds <- which(r2_plot == 0)
replace_inds2 <- which(r2_plot == 0.01)
replace_inds3 <- which(r2_plot == 0.98)
r2_plot <- as.character(r2_plot)
r2_plot[replace_inds] <- "0.000" # replace the 0 with 0.000 (as R^2 was calculated to 3dp)
r2_plot[replace_inds2] <- "0.010" # add terminal 0
r2_plot[replace_inds3] <- "0.980" # add terminal 0
r2_plot <- paste0(" = ", r2_plot)

# Make a dataframe of variables r^2 values, vectors for plot placement and equation for the plot
e_eq <- data.frame(variable = eq_var_list, x_pos = rep(0.0, 24) , y1_pos = rep(1.1,24), y2_pos = rep(c(4.2, 0.1, 1.05, 0.45, 0.5, 1), each = 4), age = rep(c(0.05, 0.1, 0.5, 1.0), 6),
                   r2_x_pos = rep(0.6, 24), rsq = r2_plot, rsquare = "R^2 ",
                   group = factor(eq_var_list,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), 
                                  ordered = TRUE, 
                                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Prop. recomb. trip.","(3SEQ)")),
                                             expression(atop("Prop. res. quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) ) ) 

p <- ggplot(e, aes(x = number_of_events, y = value)) +
  geom_point(colour = "gray55", alpha = 0.2) +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
  facet_grid(group~age, scale = "free_y", labeller = label_parsed) +
  scale_x_continuous(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) + 
  geom_text(data = e_eq, aes(x = x_pos, y = y2_pos, label = rsquare), parse = TRUE, hjust = 0, size = 3) + 
  geom_text(data = e_eq, aes(x = r2_x_pos, y = y2_pos, label = rsq), parse = FALSE, hjust = 0, size = 3)
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_increasingNumEvents_allTreeDepths_regression_freey.png"), plot = p, units = "in", height = 9, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_increasingNumEvents_allTreeDepths_regression_freey.pdf"), height = 9, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = number_of_events, y = value)) +
  geom_point(colour = "gray55", alpha = 0.2) +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
  facet_grid(group~age, scale = "free_y", labeller = label_parsed) +
  scale_x_continuous(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) + 
  geom_text(data = e_eq, aes(x = x_pos, y = y2_pos, label = rsquare), parse = TRUE, hjust = 0, size = 3) + 
  geom_text(data = e_eq, aes(x = r2_x_pos, y = y2_pos, label = rsq), parse = FALSE, hjust = 0, size = 3)
dev.off()

patchwork <- p
cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_increasingNumEvents_allTreeDepths_regression_freey_withTitle.pdf"), height = 9, width = 9, fallback_resolution = 300)
patchwork + 
  plot_annotation(title = "Total tree depth (substitutions per site)",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(family = "")))
dev.off()

# Continuation of plot 4, but includes both event types
e = subset(ts2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_event_type != "none")
e = subset(e, tree_age == tree_length)
e$tree2_event_type = factor(e$tree2_event_type, ordered = TRUE)
e$type = paste(e$tree2_event_type, e$tree2_event_position)
e = e[e$variable %in% c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group <- factor(e$variable, levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Prop. recomb. trip.","(3SEQ)")),
                             expression(atop("Prop. res. quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

r2_raw <- c()
var_list <- c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed")
eq_var_list <- rep(var_list, each = 2)
event_list <- c("nonreciprocal","reciprocal")
for (var in var_list){
  var_df <- e[(e$variable == var),]
  for (ev in event_list){
    # get only entries with this variable
    # insert into regression model
    age_df <- var_df[(var_df$tree2_event_type == ev),]
    m <- lm(value ~ number_of_events, age_df)
    r2_var <- summary(m)$r.squared
    print(paste0("Variable: ",var," - event type: ",ev," - R^2 : ", round(r2_var, digits = 3) ))
    r2_raw <- c(r2_raw, r2_var) 
  }
}

# Quick function to reformat r^2 values the way I like
reformat.r2 <- function(number){
  number <- round(number, digits = 3)
}

# Apply function to r^2 values
r2_plot <- unlist(lapply(r2_raw,reformat.r2))
replace_inds <- which(r2_plot == 0)
replace_inds2 <- which(r2_plot == 0.98)
r2_plot <- as.character(r2_plot)
r2_plot[replace_inds] <- "0.000" # replace the 0 with 0.000 (as R^2 was calculated to 3dp)
r2_plot[replace_inds2] <- "0.980" # pad the 0.98 out with a 0 so it's to 3dp
r2_plot <- paste0(" = ", r2_plot)


# Make a dataframe of variables r^2 values, vectors for plot placement and equation for the plot
e_eq <- data.frame(variable = eq_var_list, x_pos = rep(0.0, 12) , y1_pos = rep(1.1,12), 
                   y2_pos = rep(c(1.67, 0.003, 1.02, 0.28, 0.28, 1.05), each = 2),
                   r2_x_pos = rep(0.5, 12), rsq = r2_plot, rsquare = "R^2 ", tree2_event_type = rep(c("nonreciprocal","reciprocal"),6),
                   fac_event_type = factor(rep(event_list,6), ordered = TRUE, labels = c("Nonreciprocal","Reciprocal")),
                   group = factor(eq_var_list,levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), 
                                  ordered = TRUE, 
                                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Prop. recomb. trip.","(3SEQ)")),
                                             expression(atop("Prop. res. quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) ) ) 
# Change e$tree2_event_type to a new factor with pretty formatting for output
e$fac_event_type <- factor(e$tree2_event_type, ordered = TRUE, labels = c("Nonreciprocal","Reciprocal"))

p <- ggplot(e, aes(x = number_of_events, y = value)) +
  geom_point(colour = "gray55", alpha = 0.2) +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
  facet_grid(group ~ fac_event_type, scale = "free_y", labeller = label_parsed) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  geom_text(data = e_eq, aes(x = r2_x_pos, y = y2_pos, label = rsq), parse = FALSE, hjust = 0, size = 3) +
  geom_text(data = e_eq, aes(x = x_pos, y = y2_pos, label = rsquare), parse = TRUE, hjust = 0, size = 3) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_increasingNumEvents_RecipNonrecip_regression_freey.png"), plot = p, units = "in", height = 7, width = 5)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_increasingNumEvents_RecipNonrecip_regression_freey.pdf"), height = 7, width = 5, fallback_resolution = 300)
ggplot(e, aes(x = number_of_events, y = value)) +
  geom_point(colour = "gray55", alpha = 0.2) +
  geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
  facet_grid(group ~ fac_event_type, scale = "free_y", labeller = label_parsed) +
  scale_x_continuous(name = "\n Proportion of DNA introgressed \n") +
  ylab("\n Test statistic value \n") +
  geom_text(data = e_eq, aes(x = r2_x_pos, y = y2_pos, label = rsq), parse = FALSE, hjust = 0, size = 3) +
  geom_text(data = e_eq, aes(x = x_pos, y = y2_pos, label = rsquare), parse = TRUE, hjust = 0, size = 3) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")))
dev.off()

patchwork <- p
cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_increasingNumEvents_RecipNonrecip_regression_freey_withTitle.pdf"), height = 7, width = 5, fallback_resolution = 300)
patchwork + 
  plot_annotation(title = "Type of introgression event",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(family = "")))
dev.off()


# Plot 5: How does reciprocity of events influence detection of treelikeness?
print("Plot 5")
e = subset(ts2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == tree_length)
e = subset(e, tree2_event_type != "none")
e = e[e$variable %in% c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group <- factor(e$variable, levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Proportion of recombinant triplets","(3SEQ)")),
                             expression(atop("Proportion of resolved quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

p <- ggplot(e, aes(x = number_of_events, y = value, colour = tree2_event_type)) +
  geom_smooth(method = "lm") +
  facet_wrap(~group, scales = "free_y", labeller = label_parsed, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n", breaks = c(0:8), labels = c(0:8)) +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 8),
        legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + 
  guides(color = guide_legend(title = "Event type")) +
  scale_colour_manual(values = c("gray60","black"),labels = c("Nonreciprocal", "Reciprocal"))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_ReciprocialAndNonreciprocalEvents_lm_freey.png"), plot = p, units = "in", height = 6.6, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_ReciprocialAndNonreciprocalEvents_lm_freey.pdf"), height = 6.6, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = number_of_events, y = value, colour = tree2_event_type)) +
  geom_smooth(method = "lm") +
  facet_wrap(~group, scales = "free_y", labeller = label_parsed, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n", breaks = c(0:8), labels = c(0:8)) +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 8),
        legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + 
  guides(color = guide_legend(title = "Event type")) +
  scale_colour_manual(values = c("gray60","black"),labels = c("Nonreciprocal", "Reciprocal"))
dev.off()

p <- ggplot(e, aes(x = number_of_events, y = value, colour = tree2_event_type)) +
  geom_smooth(method = "lm") +
  facet_wrap(~group, labeller = label_parsed, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n", breaks = c(0:8), labels = c(0:8)) +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + 
  guides(color = guide_legend(title = "Event type")) +
  scale_colour_manual(values = c("gray60","black"),labels = c("Nonreciprocal", "Reciprocal"))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_ReciprocialAndNonreciprocalEvents_lm_fixedy.png"), plot = p, units = "in", height = 6.6, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_ReciprocialAndNonreciprocalEvents_lm_fixedy.pdf"), height = 6.6, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = number_of_events, y = value, colour = tree2_event_type)) +
  geom_smooth(method = "lm") +
  facet_wrap(~group, labeller = label_parsed, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n", breaks = c(0:8), labels = c(0:8)) +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + 
  guides(color = guide_legend(title = "Event type")) +
  scale_colour_manual(values = c("gray60","black"),labels = c("Nonreciprocal", "Reciprocal"))
dev.off()

e = subset(ts3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == tree_length)
e = subset(e, tree2_event_type != "none")
e = e[e$variable %in% c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group <- factor(e$variable, levels = c("PHI_observed","proportion_recombinant_triplets","prop_resolved_quartets","mean_delta_q","mode_delta_q","neighbour_net_trimmed"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Proportion of recombinant triplets","(3SEQ)")),
                             expression(atop("Proportion of resolved quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

p <- ggplot(e, aes(x = proportion_tree2, y = value, colour = tree2_event_type)) +
  geom_smooth(method = "lm") +
  facet_wrap(~group, scales = "free_y", labeller = label_parsed, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 8),
        legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + 
  guides(color = guide_legend(title = "Event type")) +
  scale_colour_manual(values = c("gray60","black"),labels = c("Nonreciprocal", "Reciprocal"))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_ReciprocialAndNonreciprocalEvents_lm_freey.png"), plot = p, units = "in", height = 6.6, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_ReciprocialAndNonreciprocalEvents_lm_freey.pdf"), height = 6.6, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = proportion_tree2, y = value, colour = tree2_event_type)) +
  geom_smooth(method = "lm") +
  facet_wrap(~group, scales = "free_y", labeller = label_parsed, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 8),
        legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + 
  guides(color = guide_legend(title = "Event type")) +
  scale_colour_manual(values = c("gray60","black"),labels = c("Nonreciprocal", "Reciprocal"))
dev.off()

p <- ggplot(e, aes(x = proportion_tree2, y = value, colour = tree2_event_type)) +
  geom_smooth(method = "lm") +
  facet_wrap(~group, labeller = label_parsed, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + 
  guides(color = guide_legend(title = "Event type")) +
  scale_colour_manual(values = c("gray60","black"),labels = c("Nonreciprocal", "Reciprocal"))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_ReciprocialAndNonreciprocalEvents_lm_fixedy.png"), plot = p, units = "in", height = 6.6, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_ReciprocialAndNonreciprocalEvents_lm_fixedy.pdf"), height = 6.6, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = proportion_tree2, y = value, colour = tree2_event_type)) +
  geom_smooth(method = "lm") +
  facet_wrap(~group, labeller = label_parsed, ncol = 3) +
  scale_x_continuous(name = "\n Number of introgression events \n") +
  ylab("\n Test statistic value \n") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + 
  guides(color = guide_legend(title = "Event type")) +
  scale_colour_manual(values = c("gray60","black"),labels = c("Nonreciprocal", "Reciprocal"))
dev.off()




# Plot 6: increasing proportion of introgressed DNA/increasing number of events - are the results statistically significant?
# exp 3: histograms for p-values, for 0% - 50% introgression in 10% increments, for all test statistics
print("Plot 6")
e = subset(bs3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == tree_length)
e = subset(e, tree2_event_type != "none")
e = subset(e, tree2_event_type != "reciprocal")
e = e[e$variable %in% c("PHI_p_value","X3Seq_p_value","likelihood_mapping_p_value","mean_delta_q_p_value","mode_delta_q_p_value","neighbour_net_trimmed_p_value"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group <- factor(e$variable, levels = c("PHI_p_value","X3Seq_p_value","likelihood_mapping_p_value","mean_delta_q_p_value","mode_delta_q_p_value","neighbour_net_trimmed_p_value"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Prop. recomb. trip.","(3SEQ)")),
                             expression(atop("Prop. res. quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

p <- ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~proportion_tree2,scales = "free_y", labeller = label_parsed) +
  xlab("\n P value \n") +
  ylab("\n Count \n") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        strip.text.y = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1) )
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_StatisticalSignificance_freey.png"), plot = p, units = "in", height = 9, width = 9)


cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_StatisticalSignificance_freey.pdf"), height = 9, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~proportion_tree2,scales = "free_y", labeller = label_parsed) +
  xlab("\n P value \n") +
  ylab("\n Count \n") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        strip.text.y = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1) )
dev.off()

p <- ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~proportion_tree2, labeller = label_parsed) +
  xlab("\n P value \n") +
  ylab("\n Count \n") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 7), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 7), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        strip.text.y = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1) )
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_StatisticalSignificance_fixedy.png"), plot = p, units = "in", height = 9, width = 9)

patchwork <- p
cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_StatisticalSignificance_fixedy_withTitle.pdf"), height = 9, width = 9, fallback_resolution = 300)
patchwork + 
  plot_annotation(title = "Proportion of introgressed DNA (%)",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(family = "")))
dev.off()

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_StatisticalSignificance_fixedy.pdf"), height = 9, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~proportion_tree2, labeller = label_parsed) +
  xlab("\n P value \n") +
  ylab("\n Count \n") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 7), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 7), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        strip.text.y = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1) )
dev.off()

# exp 2: histograms for p-values, for 0-8 events increasing by 1 event each time, for all test statistics
e = subset(bs2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == tree_length)
e = subset(e, tree2_event_type != "reciprocal")
e = e[e$variable %in% c("PHI_p_value","X3Seq_p_value","likelihood_mapping_p_value","mean_delta_q_p_value","mode_delta_q_p_value","neighbour_net_trimmed_p_value"),]
# Have to reorder variables so the grid comes out in the right way - do this using a new column that's a factor
e$group <- factor(e$variable, levels = c("PHI_p_value","X3Seq_p_value","likelihood_mapping_p_value","mean_delta_q_p_value","mode_delta_q_p_value","neighbour_net_trimmed_p_value"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Prop. recomb. trip.","(3SEQ)")),
                             expression(atop("Prop. res. quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

p <- ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~number_of_events,scales = "free_y", labeller = label_parsed) +
  xlab("\n P value \n") +
  ylab("\n Count \n") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        strip.text.y = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_StatisticalSignificance_freey.png"), plot = p, units = "in", height = 9, width = 9)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_StatisticalSignificance_freey.pdf"), height = 9, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~number_of_events,scales = "free_y", labeller = label_parsed) +
  xlab("\n P value \n") +
  ylab("\n Count \n") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        strip.text.y = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1))
dev.off()

p <- ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~number_of_events, labeller = label_parsed) +
  xlab("\n P value \n") +
  ylab("\n Count \n") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 7), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 7), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        strip.text.y = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1))
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_StatisticalSignificance_fixedy.png"), plot = p, units = "in", height = 9, width = 9)

patchwork <- p
cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_StatisticalSignificance_fixedy_withTitle.pdf"), height = 9, width = 9, fallback_resolution = 300)
patchwork + 
  plot_annotation(title = "Number of introgression events",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(family = "")))
dev.off()


cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_StatisticalSignificance_fixedy.pdf"), height = 9, width = 9, fallback_resolution = 300)
ggplot(e, aes(x = value)) +
  geom_histogram(breaks = c(seq(0,1,0.05))) +
  facet_grid(group~number_of_events, labeller = label_parsed) +
  xlab("\n P value \n") +
  ylab("\n Count \n") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 7), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 7), strip.text = element_text(size = 10), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        strip.text.y = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  scale_x_continuous(labels = c(0,0.25,0.5,0.75,1))
dev.off()




# Plot 7: increasing proportion of introgressed DNA - are the results statistically significant?
# for exp 3: Line graph for p-values, for 0% - 50% introgression in 10% increments, for all test statistics
print("Plot 7")
e = subset(bs3_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == tree_length)
e = subset(e, tree2_event_type != "reciprocal")
e = subset(e, tree2_event_type != "none")
e = subset(e, variable != "PHI_observed_p_value")
e = subset(e, variable != "num_recombinant_sequences_p_value")
e = e[e$variable %in% c("PHI_p_value","X3Seq_p_value","likelihood_mapping_p_value","mean_delta_q_p_value","mode_delta_q_p_value","neighbour_net_trimmed_p_value"),]

# to make new df for the plot
PHI_p_value <- c(nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                 nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
X3SEQ_p_value <- c(nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                   nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
likelihood_mapping_p_value <- c(nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                                nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
mean_delta_q_p_value <- c(nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                          nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
mode_delta_q_p_value <- c(nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                          nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))
neighbour_net_trimmed_p_value <- c(nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.1,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.2,]),
                                   nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.3,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.4,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$proportion_tree2 == 0.5,]))


value <- c(PHI_p_value, X3SEQ_p_value, likelihood_mapping_p_value, mean_delta_q_p_value, mode_delta_q_p_value, neighbour_net_trimmed_p_value)
ts <- c(rep("PHI_p_value",6), rep("X3SEQ_p_value",6), rep("likelihood_mapping_p_value",6), rep("mean_delta_q_p_value",6), rep("mode_delta_q_p_value",6), rep("neighbour_net_trimmed_p_value",6))
proportion_introgressed_DNA <- c(rep(seq(0,0.5,0.1),6))
f <- data.frame(proportion_introgressed_DNA,ts,value,stringsAsFactors = FALSE)
f$variable <- factor(ts, levels = c("PHI_p_value","X3SEQ_p_value","likelihood_mapping_p_value","mean_delta_q_p_value","mode_delta_q_p_value","neighbour_net_trimmed_p_value"), ordered = TRUE, 
                  labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Proportion of recombinant triplets","(3SEQ)")),
                             expression(atop("Proportion of resolved quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                             expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

p <- ggplot(f, aes(x = proportion_introgressed_DNA, y = value)) +
  geom_line(size=0.5) +
  facet_wrap(~variable, labeller = label_parsed, nrow = 3, ncol = 3) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 8),
        strip.text.y = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 10)) +
  scale_x_continuous(name = "\n Proportion of introgressed DNA (%) \n", labels = seq(0,0.5,0.1), breaks = seq(0,0.5,0.1), minor_breaks = c(), limits = c(0,0.5)) + 
  scale_y_continuous(name = "\n Percent of simulations that reject the null hypothesis (%) \n (p-value < 0.05) \n",
                     labels = seq(0,100,10), breaks = seq(0,100,10), minor_breaks = seq(0,100,5), limits = c(0,100)) + 
  geom_hline(aes(yintercept = 5, colour = "red"), linetype = "dashed", size = 0.5) + 
  scale_colour_manual("Ideal false\npositive rate\n", values="red", labels = "5% when\n\u03b1 = 0.05")
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_statisticalSignificance_summaryLines_fixedy.png"), plot = p, units = "in", width = 7.5, height = 6)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp3_statisticalSignificance_summaryLines_fixedy.pdf"), width = 7.5, height = 6, fallback_resolution = 300)
ggplot(f, aes(x = proportion_introgressed_DNA, y = value)) +
  geom_line(size=0.5) +
  facet_wrap(~variable, labeller = label_parsed, nrow = 3, ncol = 3) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 8),
        strip.text.y = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 10)) +
  scale_x_continuous(name = "\n Proportion of introgressed DNA (%) \n", labels = seq(0,0.5,0.1), breaks = seq(0,0.5,0.1), minor_breaks = c(), limits = c(0,0.5)) + 
  scale_y_continuous(name = "\n Percent of simulations that reject the null hypothesis (%) \n (p-value < 0.05) \n",
                     labels = seq(0,100,10), breaks = seq(0,100,10), minor_breaks = seq(0,100,5), limits = c(0,100)) + 
  geom_hline(aes(yintercept = 5, colour = "red"), linetype = "dashed", size = 0.5) + 
  scale_colour_manual("Ideal false\npositive rate\n", values="red", labels = "5% when\n\u03b1 = 0.05")
dev.off()

# and for exp2
e = subset(bs2_df, tree1_tree_shape == 'balanced')
e = subset(e, tree2_tree_shape == 'balanced')
e = subset(e, tree_age == tree_length)
e = subset(e, tree2_event_type != "reciprocal")
e = subset(e, tree2_event_type != "none")
e = subset(e, variable != "PHI_observed_p_value")
e = subset(e, variable != "num_recombinant_sequences_p_value")
e = e[e$variable %in% c("PHI_p_value","X3Seq_p_value","likelihood_mapping_p_value","mean_delta_q_p_value","mode_delta_q_p_value","neighbour_net_trimmed_p_value"),]

e_none = subset(bs2_df, tree1_tree_shape == 'balanced')
e_none = subset(e_none, tree2_tree_shape == 'balanced')
e_none = subset(e_none, tree_age == tree_length)
e_none = subset(e_none, tree2_event_type == "none")
e_none = e_none[e_none$variable %in% c("PHI_p_value","X3Seq_p_value","likelihood_mapping_p_value","mean_delta_q_p_value","mode_delta_q_p_value","neighbour_net_trimmed_p_value"),]

# to make new df for the plot
PHI_p_value <- c(nrow(e_none[e_none$variable == "PHI_p_value" & e_none$value <= 0.05 & e_none$number_of_events == 0,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 1,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 2,]),
                 nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 3,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 4,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 5,]),
                 nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 6,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 7,]),nrow(e[e$variable == "PHI_p_value" & e$value <= 0.05 & e$number_of_events == 8,]))
X3SEQ_p_value <- c(nrow(e_none[e_none$variable == "X3Seq_p_value" & e_none$value <= 0.05 & e_none$number_of_events == 0,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 1,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 2,]),
                   nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 3,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 4,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 5,]),
                   nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 6,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 7,]),nrow(e[e$variable == "X3Seq_p_value" & e$value <= 0.05 & e$number_of_events == 8,]))
likelihood_mapping_p_value <- c(nrow(e_none[e_none$variable == "likelihood_mapping_p_value" & e_none$value <= 0.05 & e_none$number_of_events == 0,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 1,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 2,]),
                                nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 3,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 4,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 5,]),
                                nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 6,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 7,]),nrow(e[e$variable == "likelihood_mapping_p_value" & e$value <= 0.05 & e$number_of_events == 8,]))
mean_delta_q_p_value <- c(nrow(e_none[e_none$variable == "mean_delta_q_p_value" & e_none$value <= 0.05 & e_none$number_of_events == 0,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 1,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 2,]),
                          nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 3,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 4,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 5,]),
                          nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 6,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 7,]),nrow(e[e$variable == "mean_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 8,]))
mode_delta_q_p_value <- c(nrow(e_none[e_none$variable == "mode_delta_q_p_value" & e_none$value <= 0.05 & e_none$number_of_events == 0,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 1,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 2,]),
                          nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 3,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 4,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 5,]),
                          nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 6,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 7,]),nrow(e[e$variable == "mode_delta_q_p_value" & e$value <= 0.05 & e$number_of_events == 8,]))
neighbour_net_trimmed_p_value <- c(nrow(e_none[e_none$variable == "neighbour_net_trimmed_p_value" & e_none$value <= 0.05 & e_none$number_of_events == 0,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 1,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 2,]),
                                   nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 3,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 4,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 5,]),
                                   nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 6,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 7,]),nrow(e[e$variable == "neighbour_net_trimmed_p_value" & e$value <= 0.05 & e$number_of_events == 8,]))

value <- c(PHI_p_value, X3SEQ_p_value, likelihood_mapping_p_value, mean_delta_q_p_value, mode_delta_q_p_value, neighbour_net_trimmed_p_value)
value = value/20*100 # 20 replicates performed for each event - divide by 20 to get percentage
ts <- c(rep("PHI_p_value",9), rep("X3SEQ_p_value",9), rep("likelihood_mapping_p_value",9), rep("mean_delta_q_p_value",9), rep("mode_delta_q_p_value",9), rep("neighbour_net_trimmed_p_value",9))
num_events <- c(rep(seq(0,8,1),6))
f <- data.frame(num_events,ts,value,stringsAsFactors = FALSE)
f$variable <- factor(ts, levels = c("PHI_p_value","X3SEQ_p_value","likelihood_mapping_p_value","mean_delta_q_p_value","mode_delta_q_p_value","neighbour_net_trimmed_p_value"), ordered = TRUE, 
                     labels = c(expression(atop("PHI","(PhiPack)")), expression(atop("Proportion of recombinant triplets","(3SEQ)")),
                                expression(atop("Proportion of resolved quartets","(IQ-Tree)")), expression(atop(paste('Mean ', delta["q"]),paste("(", delta," plots)"))),
                                expression(atop(paste('Mode ', delta["q"]),paste("(", delta," plots)"))), expression(atop("Tree proportion","(this paper)")) ) )

p <- ggplot(f, aes(x = num_events, y = value)) +
  geom_line(size=0.5) +
  facet_wrap(~variable, labeller = label_parsed, nrow = 3, ncol = 3) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 8),
        strip.text.y = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  scale_x_continuous(name = "\n Number of introgression events \n", labels = seq(0,8,1), breaks = seq(0,8,1), limits = c(0,8), minor_breaks = c()) + 
  scale_y_continuous(name = "\n Percent of simulations that reject the null hypothesis (%) \n (p-value < 0.05) \n",
                     labels = seq(0,100,10), breaks = seq(0,100,10), minor_breaks = seq(0,100,5), limits = c(0,100)) + 
  geom_hline(aes(yintercept = 5, colour = "red"), linetype = "dashed", size = 0.5) + 
  scale_colour_manual("Ideal false\npositive rate\n", values="red", labels = "5% when\n\u03b1 = 0.05")
ggsave(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_statisticalSignificance_summaryLines_fixedy.png"), plot = p, units = "in", width = 7.5, height = 6)

cairo_pdf(filename = paste0(plots_folder,run_id,"_td",tree_length,"_exp2_statisticalSignificance_summaryLines_fixedy.pdf"), width = 7.5, height = 6, fallback_resolution = 300)
ggplot(f, aes(x = num_events, y = value)) +
  geom_line(size=0.5) +
  facet_wrap(~variable, labeller = label_parsed, nrow = 3, ncol = 3) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8), strip.text = element_text(size = 8), strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 8),
        strip.text.y = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  scale_x_continuous(name = "\n Number of introgression events \n", labels = seq(0,8,1), breaks = seq(0,8,1), limits = c(0,8), minor_breaks = c()) + 
  scale_y_continuous(name = "\n Percent of simulations that reject the null hypothesis (%) \n (p-value < 0.05) \n",
                     labels = seq(0,100,10), breaks = seq(0,100,10), minor_breaks = seq(0,100,5), limits = c(0,100)) + 
  geom_hline(aes(yintercept = 5, colour = "red"), linetype = "dashed", size = 0.5) + 
  scale_colour_manual("Ideal false\npositive rate\n", values="red", labels = "5% when\n\u03b1 = 0.05")
dev.off()
