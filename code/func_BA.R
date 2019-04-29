empirical.runTS <- function(alignment_path, program_paths){
  # Open the nexus file and get the number of taxa and the number of characters 
  n <- read.nexus.data(alignment_path)
  n_taxa <- length(n)
  n_char <- length(unlist(n[1]))
  
  # Run IQ-tree on the alignment (if it hasn't already been run), and get the likelihood mapping results
  call.IQTREE.quartet(program_paths[["IQTree"]],alignment_path,n_taxa)
  
  # output_dataframe: dataset, loci name, number of taxa, number of characters, test statistic values
}

empirical.parametric.bootstrap <- function(){
  # If it hasn't already been run, call and run IQTree
  call.IQTREE(program_paths["IQTree"],alignment_path)
  
  #Extract the parameters from the .iqtree log file.
  params <- get.simulation.parameters(paste0(alignment_path,".iqtree"))
  
  # output_dataframe: dataset, loci name, number of taxa, number of characters, test statistic values
}