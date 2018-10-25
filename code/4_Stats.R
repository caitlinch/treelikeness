# R code to run statistical tests on the none event/1 event data
# Specify which file paths to use

run_location = "mac"
# run_location = "soma"

if (run_location == "mac"){
  # Set file paths etc
  input_folder <- "/Users/caitlincherryh/Documents/Results/simulations_20180913/collatedOutput/"
  output_folder <- "/Users/caitlincherryh/Documents/Results/simulations_20180913/plots/plots_20181017/"
  # Set working directory
  maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # for work computer
} else if (run_location == "soma") {
  # Set file paths etc
  input_folder <- "/data/caitlin/treelikeness/results/"
  
  # Set working directory
  maindir <- "/data/caitlin/treelikeness/"
}

# load required libraries


# Open dataframes of data
plot1_df <- read.csv(paste0(input_folder,"plot1_testStatistics_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)
plot2_df <- read.csv(paste0(input_folder,"plot2_testStatistics_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)
plot3_df <- read.csv(paste0(input_folder,"plot3_testStatistics_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)
plot4_df <- read.csv(paste0(input_folder,"plot4_testStatistics_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)
bs_df <-  read.csv(paste0(input_folder,"plot4_p_value_collatedSimulationData_melted.csv"), stringsAsFactors = FALSE)

# create empty vectors for storage
test_statistic<- c()
event1 <- c()
event1_mean <- c()
event2 <- c()
event2_mean <- c()
t <- c()
df <- c()
p_value <- c()
null_value <- c()

# Get the different test statistics
vars <- unique(plot1_df$variable)
# for each test statistic, extract the list of values for each event
for (var in vars){
  print(var)
  sub_df <- subset(plot1_df,variable == var)
  # Extract the test statistic values for each event type
  none <- subset(sub_df, tree2_event_position == "none")
  none <- subset(none, tree2_event_type == "none")
  close_recip <- subset(sub_df, tree2_event_position == "close")
  close_recip <- subset(close_recip, tree2_event_type == "reciprocal")
  close_nonrecip <- subset(sub_df, tree2_event_position == "close")
  close_nonrecip <- subset(close_nonrecip, tree2_event_type == "nonreciprocal")
  diver_recip <- subset(sub_df, tree2_event_position == "divergent")
  diver_recip <- subset(diver_recip, tree2_event_type == "reciprocal")
  diver_nonrecip <- subset(sub_df, tree2_event_position == "divergent")
  diver_nonrecip <- subset(diver_nonrecip, tree2_event_type == "nonreciprocal")
  anci_recip <- subset(sub_df, tree2_event_position == "ancient")
  anci_recip <- subset(anci_recip, tree2_event_type == "reciprocal")
  anci_nonrecip <- subset(sub_df, tree2_event_position == "ancient")
  anci_nonrecip <- subset(anci_nonrecip, tree2_event_type == "nonreciprocal")
  # Make a list of event type dfs
  event_dfs <- list("close_recip" = close_recip,"close_nonrecip" = close_nonrecip,"diver_recip" = diver_recip,"diver_nonrecip" = diver_nonrecip,"anci_recip" = anci_recip,"anci_nonrecip" = anci_nonrecip)
  event_names <- c("close_recip","close_nonrecip","diver_recip","diver_nonrecip","anci_recip","anci_nonrecip")
  # iterate through event types
  for (event in event_names){
    print(event)
    # conduct the t test
    t_test <- t.test(none$value, event_dfs[[event]]$value, var.equal = FALSE)
    # extract the values
    test_statistic <- c(test_statistic,var)
    event1 <- c(event1,"none")
    event1_mean <- c(event1_mean,t_test[["estimate"]][1])
    event2 <- c(event2,event)
    event2_mean <- c(event2_mean,t_test[["estimate"]][2])
    t <- c(t,t_test[["statistic"]])
    df <- c(df, t_test[["parameter"]])
    p_value <- c(p_value, t_test[["p.value"]])
    null_value <- c(null_value,t_test[["null.value"]])
  }
  print("")
  }

stats_df <- data.frame(test_statistic,event1,event1_mean,event2,event2_mean,t,df,p_value,null_value)
stats_df <- stats_df[order(stats_df$p_value),]
write.csv(stats_df, file = paste0(input_folder,"eventType_t_tests.csv"))

