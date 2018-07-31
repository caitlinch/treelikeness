# Quick script to open the full .csv files with the numerators and denominators fractions for test statistics 1, 4 and 5, and create some quick plots

library(reshape2)
library(ggplot2)
library(gridExtra)

output_folder <- "/Users/caitlincherryh/Documents/TestAlignmentResults/0_prelim_mk3/"
internal_df <- read.csv("/Users/caitlincherryh/Documents/TestAlignmentResults/0_prelim_mk3/internal_collated_output_includesFracs_fixedPDM.csv")
external_df <- read.csv("/Users/caitlincherryh/Documents/TestAlignmentResults/0_prelim_mk3/external_collated_output_includesFracs_fixedPDM.csv")
phylo_df <- read.csv("/Users/caitlincherryh/Documents/TestAlignmentResults/0_prelim_mk3/2trees_collated_output_includesFracs_fixedPDM.csv")
phylo_df["proportion_tree2"] <- 1-phylo_df$proportion_tree1

#external_melt <- melt(external_df,id.vars = c(4),measure.vars = c(23,24,25,26,27,28))
#internal_melt <- melt(internal_df,id.vars = c(3),measure.vars = c(23,24,25,26,27,28))
#phylo_melt <- melt(phylo_df,id.vars = c(32),measure.vars = c(26,27,28,29,30,31))

#external
external_melt1 <- melt(external_df,id.vars = c(4,23,24,25,26,27,28), measure.vars = c())
ex_ts1 <- ggplot(external_melt1, aes(x=external_recombination)) + geom_point(size = 2, col = "#009E73", shape = 16, aes(y=ts1_num)) + 
  geom_point(size = 2, col = "#E69F00", shape = 17, aes(y=ts1_denom)) + stat_summary(fun.y = mean, geom="line", lwd=0.5, col = "#009E73", aes(x=external_recombination, y=ts1_num)) +
  stat_summary(fun.y = mean, geom="line", lwd=0.5, col = "#E69F00", aes(x=external_recombination, y=ts1_denom)) + labs(x = "External recombination (%)", y = "Statistic value", title = "TS1: Splittable Percentage")

ex_ts4 <- ggplot(external_melt1, aes(x=external_recombination)) + geom_point(size = 2, col = "#009E73", shape = 16, aes(y=ts4_num)) + 
  geom_point(size = 2, col = "#E69F00", shape = 17, aes(y=ts4_denom)) + stat_summary(fun.y = mean, geom="line", lwd=0.5, col = "#009E73", aes(x=external_recombination, y=ts4_num)) +
  stat_summary(fun.y = mean, geom="line", lwd=0.5, col = "#E69F00", aes(x=external_recombination, y=ts4_denom)) + labs(x = "External recombination (%)", y = "Statistic value", title = "TS4: Split Decomposition")

ex_ts5 <- ggplot(external_melt1, aes(x=external_recombination)) + geom_point(size = 2, shape = 16, aes(y=ts5_num, colour = "Numerator")) + 
  geom_point(size = 2, shape = 17, aes(y=ts5_denom, colour = "Denominator")) + stat_summary(fun.y = mean, geom="line", lwd=0.5, aes(x=external_recombination, y=ts5_num,colour = "Numerator")) +
  stat_summary(fun.y = mean, geom="line", lwd=0.5, aes(x=external_recombination, y=ts5_denom, colour = "Denominator")) + 
  labs(x = "External recombination (%)", y = "Statistic value", title = "TS5: Neighbor-Net") + 
  scale_colour_manual(name="Color",values=c(Numerator = "#009E73", Denominator ="#E69F00"))

combined_plot <- grid.arrange(ex_ts1, ex_ts4, ex_ts5 , nrow = 1, ncol = 3)
plot_title <- paste0(output_folder, "externalRecombination_FractionPlots.pdf")
ggsave(filename = plot_title, plot = combined_plot, dev = "pdf", width = 32.44, height = 14.88, units = "in")
plot_title <- paste0(output_folder, "externalRecombination_FractionPlots.png")
ggsave(filename = plot_title, plot = combined_plot, dev = "png", width = 32.44, height = 14.88, units = "in")

# internal
internal_melt1 <- melt(internal_df,id.vars = c(3,23,24,25,26,27,28), measure.vars = c())
in_ts1 <- ggplot(internal_melt1, aes(x=internal_recombination)) + geom_point(size = 2, col = "#009E73", shape = 16, aes(y=ts1_num)) + 
  geom_point(size = 2, col = "#E69F00", shape = 17, aes(y=ts1_denom)) + stat_summary(fun.y = mean, geom="line", lwd=0.5, col = "#009E73", aes(x=internal_recombination, y=ts1_num)) +
  stat_summary(fun.y = mean, geom="line", lwd=0.5, col = "#E69F00", aes(x=internal_recombination, y=ts1_denom)) + labs(x = "Internal recombination (%)", y = "Statistic value", title = "TS1: Splittable Percentage")

in_ts4 <- ggplot(internal_melt1, aes(x=internal_recombination)) + geom_point(size = 2, col = "#009E73", shape = 16, aes(y=ts4_num)) + 
  geom_point(size = 2, col = "#E69F00", shape = 17, aes(y=ts4_denom)) + stat_summary(fun.y = mean, geom="line", lwd=0.5, col = "#009E73", aes(x=internal_recombination, y=ts4_num)) +
  stat_summary(fun.y = mean, geom="line", lwd=0.5, col = "#E69F00", aes(x=internal_recombination, y=ts4_denom)) + labs(x = "Internal recombination (%)", y = "Statistic value", title = "TS4: Split Decomposition")

in_ts5 <- ggplot(internal_melt1, aes(x=internal_recombination)) + geom_point(size = 2, shape = 16, aes(y=ts5_num, colour = "Numerator")) + 
  geom_point(size = 2, shape = 17, aes(y=ts5_denom, colour = "Denominator")) + stat_summary(fun.y = mean, geom="line", lwd=0.5, aes(x=internal_recombination, y=ts5_num,colour = "Numerator")) +
  stat_summary(fun.y = mean, geom="line", lwd=0.5, aes(x=internal_recombination, y=ts5_denom, colour = "Denominator")) + 
  labs(x = "Internal recombination (%)", y = "Statistic value", title = "TS5: Neighbor-Net") + 
  scale_colour_manual(name="Color",values=c(Numerator = "#009E73", Denominator ="#E69F00"))

combined_plot <- grid.arrange(in_ts1, in_ts4, in_ts5 , nrow = 1, ncol = 3)
plot_title <- paste0(output_folder, "internalRecombination_FractionPlots.pdf")
ggsave(filename = plot_title, plot = combined_plot, dev = "pdf", width = 32.44, height = 14.88, units = "in")
plot_title <- paste0(output_folder, "internalRecombination_FractionPlots.png")
ggsave(filename = plot_title, plot = combined_plot, dev = "png", width = 32.44, height = 14.88, units = "in")

# phylo
phylo_melt1 <- melt(phylo_df,id.vars = c(32,26,27,28,29,30,31), measure.vars = c())
phy_ts1 <- ggplot(phylo_melt1, aes(x=proportion_tree2)) + geom_point(size = 2, col = "#009E73", shape = 16, aes(y=ts1_num)) + 
  geom_point(size = 2, col = "#E69F00", shape = 17, aes(y=ts1_denom)) + stat_summary(fun.y = mean, geom="line", lwd=0.5, col = "#009E73", aes(x=proportion_tree2, y=ts1_num)) +
  stat_summary(fun.y = mean, geom="line", lwd=0.5, col = "#E69F00", aes(x=proportion_tree2, y=ts1_denom)) + labs(x = "Proportion of Tree 2 (%)", y = "Statistic value", title = "TS1: Splittable Percentage")

phy_ts4 <- ggplot(phylo_melt1, aes(x=proportion_tree2)) + geom_point(size = 2, col = "#009E73", shape = 16, aes(y=ts4_num)) + 
  geom_point(size = 2, col = "#E69F00", shape = 17, aes(y=ts4_denom)) + stat_summary(fun.y = mean, geom="line", lwd=0.5, col = "#009E73", aes(x=proportion_tree2, y=ts4_num)) +
  stat_summary(fun.y = mean, geom="line", lwd=0.5, col = "#E69F00", aes(x=proportion_tree2, y=ts4_denom)) + labs(x = "Proportion of Tree 2 (%)", y = "Statistic value", title = "TS4: Split Decomposition")

phy_ts5 <- ggplot(phylo_melt1, aes(x=proportion_tree2)) + geom_point(size = 2, shape = 16, aes(y=ts5_num, colour = "Numerator")) + 
  geom_point(size = 2, shape = 17, aes(y=ts5_denom, colour = "Denominator")) + stat_summary(fun.y = mean, geom="line", lwd=0.5, aes(x=proportion_tree2, y=ts5_num,colour = "Numerator")) +
  stat_summary(fun.y = mean, geom="line", lwd=0.5, aes(x=proportion_tree2, y=ts5_denom, colour = "Denominator")) + 
  labs(x = "Proportion of Tree 2 (%)", y = "Statistic value", title = "TS5: Neighbor-Net") + 
  scale_colour_manual(name="Color",values=c(Numerator = "#009E73", Denominator ="#E69F00"))

combined_plot <- grid.arrange(phy_ts1, phy_ts4, phy_ts5 , nrow = 1, ncol = 3)
plot_title <- paste0(output_folder, "phylo_FractionPlots_mean.pdf")
ggsave(filename = plot_title, plot = combined_plot, dev = "pdf", width = 32.44, height = 14.88, units = "in")
plot_title <- paste0(output_folder, "phylo_FractionPlots_mean.png")
ggsave(filename = plot_title, plot = combined_plot, dev = "png", width = 32.44, height = 14.88, units = "in")


phy_ts1 <- ggplot(phylo_melt1, aes(x=proportion_tree2)) + geom_point(size = 2, shape = 16, aes(y=ts1_num, colour = "Numerator")) + 
  geom_point(size = 2, shape = 17, aes(y=ts1_denom, colour = "Denominator")) + geom_smooth(method = "lm", aes(x = proportion_tree2, y = ts1_num, colour = "Numerator")) + 
  geom_smooth(method = "lm", aes(x = proportion_tree2, y = ts1_denom, colour = "Denominator")) + scale_colour_manual(name="Color",values=c(Numerator = "#009E73", Denominator ="#E69F00")) + 
  labs(x = "Proportion of Tree 2 (%)", y = "Statistic value", title = "TS1: Splittable Percentage")

phy_ts4 <- ggplot(phylo_melt1, aes(x=proportion_tree2)) + geom_point(size = 2, col = "#009E73", shape = 16, aes(y=ts4_num)) + 
  geom_point(size = 2, col = "#E69F00", shape = 17, aes(y=ts4_denom)) + geom_smooth(method = "lm", aes(x = proportion_tree2, y = ts4_num, colour = "Numerator")) + 
  geom_smooth(method = "lm", aes(x = proportion_tree2, y = ts4_denom, colour = "Denominator")) + scale_colour_manual(name="Color",values=c(Numerator = "#009E73", Denominator ="#E69F00"))+ 
  labs(x = "Proportion of Tree 2 (%)", y = "Statistic value", title = "TS4: Split Decomposition")

phy_ts5 <- ggplot(phylo_melt1, aes(x=proportion_tree2)) + geom_point(size = 2, shape = 16, aes(y=ts5_num, colour = "Numerator")) + 
  geom_point(size = 2, shape = 17, aes(y=ts5_denom, colour = "Denominator")) + geom_smooth(method = "lm", aes(x = proportion_tree2, y = ts5_num, colour = "Numerator")) + 
  geom_smooth(method = "lm", aes(x = proportion_tree2, y = ts5_denom, colour = "Denominator")) +
  labs(x = "Proportion of Tree 2 (%)", y = "Statistic value", title = "TS5: Neighbor-Net") + 
  scale_colour_manual(name="Color",values=c(Numerator = "#009E73", Denominator ="#E69F00"))

combined_plot <- grid.arrange(phy_ts1, phy_ts4, phy_ts5 , nrow = 1, ncol = 3)
plot_title <- paste0(output_folder, "phylo_FractionPlots_lm.pdf")
ggsave(filename = plot_title, plot = combined_plot, dev = "pdf", width = 32.44, height = 14.88, units = "in")
plot_title <- paste0(output_folder, "phylo_FractionPlots_lm.png")
ggsave(filename = plot_title, plot = combined_plot, dev = "png", width = 32.44, height = 14.88, units = "in")
