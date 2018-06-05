# Functions to create alignments for simulations with a coalescent or phylogenetic approach

library(TreeSim)
library(phytools)
library(seqinr)
library(ape)
library(phangorn)

# Create a function to make SimBac alignments
SimBac.make1 <- function(simbac_path, output_folder, ntaxa, nsites, gap, mutation_rate = 0.01, internal_recombination, external_recombination, id = ""){
  # note - site specific mutation rate defaults to 0.01 when not specified in SimBac
  # Good default gap size is 1,000,000 (1000000)
  # Create a filename
  output_file <- paste0(output_folder,"SimBac_",ntaxa,"_",nsites,"_",internal_recombination,"_",external_recombination,"_",mutation_rate,"_NA_NA_NA_",id,".fasta")
  tree_file <- paste0(output_folder,"SimBac_",ntaxa,"_",nsites,"_",internal_recombination,"_",external_recombination,"_",mutation_rate,"_NA_NA_NA_tree_",id,".treefile")
  clonal_file <- paste0(output_folder,"SimBac_",ntaxa,"_",nsites,"_",internal_recombination,"_",external_recombination,"_",mutation_rate,"_NA_NA_NA_clonalGenealogy_",id,".txt")
  arg_file <- paste0(output_folder,"SimBac_",ntaxa,"_",nsites,"_",internal_recombination,"_",external_recombination,"_",mutation_rate,"_NA_NA_NA_ancestralRecombinationGraph_",id,".gv")
  # Put together the SimBac command
  simbac_command <- paste0(simbac_path," -N ",ntaxa," -B ",nsites," -G ",gap," -T ",mutation_rate," -R ",internal_recombination,
                           " -r ",external_recombination," -o ",output_file," -c ",clonal_file," -l ",tree_file," -d ",arg_file)
  # Call SimBac
  system(simbac_command)
  
  # output a text file with all the parameters
  output_name_template <- paste0(output_folder,"SimBac_",ntaxa,"_",nsites,"_",internal_recombination,"_",external_recombination,"_",mutation_rate,"_NA_NA_NA_",id,"_params.csv") # create a name for the output file 
  row <- c("coalescent",ntaxa,nsites,internal_recombination,external_recombination,mutation_rate,"NA","NA","NA","NA","NA","NA","NA",id) # gather up all the variables
  names <- c("method","n_taxa","n_sites","internal_recombination","external_recombination","mutation_rate","birth_rate","death_rate","tree_age","mean_molecular_rate",
             "sd_molecular_rate","proportion_tree1","proportion_tree2","id") # gather up the var names
  df <- data.frame(matrix(nrow=0,ncol=14)) # make an empty dataframe
  df <- rbind(df,row) # attach the info to the empty df
  names(df) <- names # rename it so it's pretty and also actually helpful
  write.csv(df, file = output_name_template) # write the csv so you can use it later. 
}

# Function to take a row from a dataframe and separate it into its components, then call the function to make 1 SimBac alignment (and its associated tree/parameter files)
SimBac.wrapper <- function(row,program_paths){
  # Extract the path to the SimBac executable from the program_paths vector
  simbac_path <- program_paths[["SimBac"]]
  # Extract values for creating the phylogenetic alignment from the input row
  ntaxa <- as.numeric(row$n_taxa)
  nsites <- as.numeric(row$n_sites)
  gap <- as.numeric(row$gap)
  mutation_rate <- as.numeric(row$mutation_rate)
  internal_recombination <- as.numeric(row$internal_recombination)
  external_recombination <- as.numeric(row$external_recombination)
  id <- paste0(row$id,"_",row$rep)
  # Create an output folder name using 
  output_folder <- paste0(row$output_folder,"SimBac_",ntaxa,"_",nsites,"_",internal_recombination,"_",external_recombination,"_",mutation_rate,"_NA_NA_NA_",id,"/")
  # Call to phylo.make1 function to create one (1) simulation and store all information about that simulation in the folder from above
  SimBac.make1(simbac_path, output_folder, ntaxa, nsites, gap, mutation_rate, internal_recombination, external_recombination, id)
  # return the output folder so that you can open it and run the test statistics
  return(output_folder)
}

# Function to run one entire simulation using a coalescent framework : create the alignment, run test statistics, and save
SimBac.run1sim <- function(row, program_paths){
  # Call the wrapper function to create the alignment: output folder is returned so you know where to look to find the files
  sim_folder <- SimBac.wrapper(row)
}

# Create a function to make phylogenetic alignments (as outlined in simulation scheme)
# K is the proportion of the SECOND tree that will be included (provide a single value)
phylo.make1 <- function(output_folder, ntaxa, nsites, birth_rate = 0.5, tree_age = 1, mol_rate, mol_rate_sd = 0.1, K = 0,id){
  # Randomly select a death rate using a uniform distribution from 0 to 99% of the birth rate
  death_rate = runif(1,min = 0, max = (0.99*birth_rate))
  # 1. Simulate a tree
  # simulate a birth-death tree on a fixed number of extant taxa
  # n = number of taxa, numbsim = # of trees to simulate, lambda = speciation rate [good default = 0.5 from Duchenne and Lanfear (2015)]
  # mu = extinction rate (default for this project = 0)
  tree_sim <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = birth_rate, mu = death_rate, frac = 1, age = tree_age, mrca = FALSE)[[1]]
  tree_sim$edge.length <- tree_sim$edge.length * (tree_age / max(branching.times(tree_sim))) # adjust branches to exactly equal tree age (adjusts for rounding errors)
  # 2. Simulate rate variation
  # Default for mol_rate and mol_rate_sd = 0.1 as in Duchenne and Lanfear (2015)
  phylo_sim <- tree_sim
  phylo_sim$edge.length <- tree_sim$edge.length * rlnorm(length(tree_sim$edge.length), meanlog = log(mol_rate), sdlog = mol_rate_sd) # adjust branch lengths - mol rate will control the tree depth in substitutions per site
  # scale tree to have a total depth of tree age
  phylo_sim <- rescale(phylo_sim,"depth",tree_age)
  # Save the first tree and a picture of the first tree
  pdf(file = paste0(output_folder,"Phylo_",ntaxa,"_",nsites,"_NA_NA_",tree_age,"_tree1_",id,".pdf"))
  plot.phylo(phylo_sim)
  dev.off()
  write.tree(phylo_sim, file = paste0(output_folder,"Phylo_",ntaxa,"_",nsites,"_NA_NA_",tree_age,"_tree1_",id,".treefile"), tree.names = TRUE)
  
  # Perform a single SPR move to get a new tree 
  phylo_sim_2 <- rSPR(phylo_sim, moves=1) # perform a single SPR move at random
  # to get SPR distance between the trees: SPR.dist()
  # spr_dist <- sprdist
  # Save the second tree and a picture of the second tree
  pdf(file = paste0(output_folder,"Phylo_",ntaxa,"_",nsites,"_NA_NA_",tree_age,"_tree2_",id,".pdf"))
  plot.phylo(phylo_sim_2)
  dev.off()
  write.tree(phylo_sim_2, file = paste0(output_folder,"Phylo_",ntaxa,"_",nsites,"_NA_NA_",tree_age,"_tree2_",id,".treefile"), tree.names = TRUE)
  
  # Calculate how many sites of each tree will be needed 
  J <- 1 - K # proportion of first tree that will be included
  J_sites <- nsites * J # find number of sites to model from first tree
  K_sites <- nsites * K # find number of sites to model from second tree
  if ((J_sites+K_sites) < nsites){ 
    # if there are less sites then necessary, randomly add the missing sites to either tree one or tree 2
    add <- nsites - (J_sites+K_sites)
    rand <- sample(c("J","K"),1)
    if (rand == "K"){
      K_sites <- K_sites + add
    } else if (rand == "J"){
      J_sites <- J_sites + add
    }
  } else if ((J_sites+K_sites) > nsites){
    # if there are more sites then necessary, randomly subtract the missing sites from either tree one or tree 2
    subtract <- (J_sites+K_sites) - nsites
    rand <- sample(c("J","K"),1)
    if (rand == "K"){
      K_sites <- K_sites - subtract
    } else if (rand == "J"){
      J_sites <- J_sites - subtract
    }
  }
  
  # simulate the DNA along the trees, using the J and K proportions.
  dna_sim_1 <- simSeq(phylo_sim,l = K_sites) # simulate sites along the first tree
  dna_sim_2 <- simSeq(phylo_sim_2,l = J_sites) # simulate sites along the second tree
  dna_sim <- c(dna_sim_1,dna_sim_2) # concatenate the two alignments
  output_name_template <- paste0(output_folder,"Phylo_",ntaxa,"_",nsites,"_NA_NA_",tree_age,"_",K,"_",id,".nexus") # create a name for the output file
  write.phyDat(dna_sim,file = output_name_template, format = "nexus",interleaved = TRUE, datablock = FALSE) # write the output as a nexus file
  # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
  nexus <- readLines(output_name_template) # open the new nexus file
  ind <- grep("BEGIN CHARACTERS",nexus)+2 # find which line
  nexus[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
  writeLines(nexus,output_name_template) # output the edited nexus file
  
  # output a text file with all the parameters
  output_name_template <- paste0(output_folder,"Phylo_",ntaxa,"_",nsites,"_NA_NA_NA_",tree_age,"_",mol_rate,"_",K,"_",id,"_params.csv") # create a name for the output file 
  row <- c("phylogenetic",ntaxa,nsites,"NA","NA","NA",birth_rate,death_rate,tree_age,mol_rate,mol_rate_sd,J,K,id) # gather up all the variables
  names <- c("method","n_taxa","n_sites","internal_recombination","external_recombination","mutation_rate","birth_rate","death_rate","tree_age","mean_molecular_rate",
             "sd_molecular_rate","proportion_tree1","proportion_tree2","id") # gather up the var names
  df <- data.frame(matrix(nrow=0,ncol=14)) # make an empty dataframe
  df <- rbind(df,row) # attach the info to the empty df
  names(df) <- names # rename it so it's pretty and also actually helpful
  write.csv(df, file = output_name_template) # write the csv so you can use it later. 
}

# Function to take a row from a dataframe and separate it into its components, then call the function to make 1 phylogenetic alignment (and its associated tree/parameter files)
phylo.wrapper <- function(row){
  # Extract values for creating the phylogenetic alignment from the input row (convert to numeric so can use the elements for ~ maths things ~)
  ntaxa <- as.numeric(row$n_taxa)
  nsites <- as.numeric(row$n_sites)
  birth_rate <- as.numeric(row$birth_rate)
  tree_age <- as.numeric(row$tree_age)
  mol_rate <- as.numeric(row$mean_molecular_rate)
  mol_rate_sd <- as.numeric(row$sd_molecular_rate)
  K <- as.numeric(row$proportion_tree2)
  id <- paste0(row$id,"_",row$rep)
  # Create an output folder name using 
  output_folder <- paste0(row$output_folder,"Phylo_",ntaxa,"_",nsites,"_NA_NA_NA_",tree_age,"_",mol_rate,"_",K,"_",id,"/")
  # Call to phylo.make1 function to create one (1) simulation and store all information about that simulation in the folder from above
  phylo.make1(output_folder, ntaxa, nsites, birth_rate, tree_age, mol_rate, mol_rate_sd, K, id)
  # return the output folder so that you can open it and run the test statistics
  return(output_folder)
}

# Function to run one entire simulation using a phylogenetic framework : create the alignment, run test statistics, and save 
phylo.run1sim <- function(row, program_paths){
  # Call the wrapper function to create the alignment: output folder is returned so you know where to look to find the files
  sim_folder <- phylo.wrapper(row)
}











