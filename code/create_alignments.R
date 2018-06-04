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
  output_file <- paste0(output_folder,"SimBac_",ntaxa,"_",nsites,"_",internal_recombination,"_",external_recombination,"_NA_NA_",id,".fasta")
  simbac_command <- paste0(simbac_path," -N ",ntaxa," -B ",nsites," -G ",gap," -T ",mutation_rate," -R ",internal_recombination,
                           " -r ",external_recombination," -o ",output_file)
  system(simbac_command)
}

# Create a function to make phylogenetic alignments (as outlined in simulation scheme)
# K is the proportion of the SECOND tree that will be included (provide a single value)
phylo.make1 <- function(output_folder, ntaxa, nsites, birth_rate = 0.5, tree_age = 1, mol_rate = 0.1, mol_rate_sd = 0.1, K = 0,id=""){
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
  # Save the tree
  pdf(file = paste0(output_folder,"Phylo_",ntaxa,"_",nsites,"_NA_NA_",tree_age,"_tree1_",id,".pdf"))
  plot.phylo(phylo_sim)
  dev.off()
  write.tree(phylo_sim, file = paste0(output_folder,"Phylo_",ntaxa,"_",nsites,"_NA_NA_",tree_age,"_tree1_",id,".treefile"), tree.names = TRUE)
  
  # Perform a single SPR move to get a new tree 
  phylo_sim_2 <- rSPR(phylo_sim, moves=1) # perform a single SPR move at random
  
  # to get SPR distance between the trees: SPR.dist()
  spr_dist <- sprdist
  pdf(file = paste0(output_folder,"Phylo_",ntaxa,"_",nsites,"_NA_NA_",tree_age,"_tree2_",id,".pdf"))
  plot.phylo(phylo_sim_2)
  dev.off()
  write.tree(phylo_sim_2, file = paste0(output_folder,"Phylo_",ntaxa,"_",nsites,"_NA_NA_",tree_age,"_tree2_",id,".treefile"), tree.names = TRUE)
  J_vector <- 1 - K_vector # proportion of first tree that will be included
  dna_sim_1 <- as.DNAbin(simSeq(phylo_sim,l = nsites)) # simulate along the entire first tree
  dna_sim_2 <- as.DNAbin(simSeq(phylo_sim_2,l = nsites)) # simulate along the entire second tree
  output_name_template <- paste0(output_folder,"Phylo_",ntaxa,"_",nsites,"_NA_NA_",tree_age,"_")
  lapply(J_vector,mosaic.alignment,nsites,ntaxa,output_name_template,id,dna_sim_1,dna_sim_2)
}

# Want a function that, given a J value and 2 alignments (in DNAbin format), makes a concatenated alignment containing J% of tree 1 and saves it
mosaic.alignment <- function(J,nsites,ntaxa,output_name_template,id,alignment1,alignment2){
  K <- 1-J
  J_start <- 1 # find starting index for tree 1
  J_end <- floor(nsites*J) # find ending index for tree 1
  K_start <- 1 # find starting index for tree 2
  K_end <- floor(nsites*K)  # find ending index for tree 2
  if ((K_end+J_end) < nsites){
    # If there are less than 1300 base pairs due to rounding
    add <- 1300-(K_end+J_end) # work out how many base pairs to add
    rand <- sample(c("J","K"),1) # randomly pick J or K, add base pairs to one so the total is 1300
    if (rand == "K"){
      K_end <- K_end + add
    } else if (rand == "K"){
      J_end <- J_end + add
    }
  }
  if ((K_end+J_end) > nsites){
    # If there are less than 1300 base pairs due to rounding
    subtract <- (K_end+J_end)-1300 # work out how many base pairs to subtract
    rand <- sample(c("J","K"),1) # randomly pick J or K, remove base pairs from one so the total is 1300
    if (rand == "K"){
      K_end <- K_end - subtract
    } else if (rand == "K"){
      J_end <- J_end - subtract
    }
  }
  # Create a new mosaic alignment using the start and end indices
  alignment1_concat <- alignment1[1:ntaxa,J_start:J_end] # get the proportion of first alignment
  alignment2_concat <- alignment2[1:ntaxa,K_start:K_end] # get the proportion of second alignment
  dna_sim <- cbind(alignment1_concat,alignment2_concat) # concatenate the alignments
  output_name <- paste0(output_name_template,"K",K,"_",id,".nexus") # create a name for the output file 
  write.nexus.data(dna_sim,file=output_name,format = "dna",interleaved = TRUE, datablock = FALSE) # write the output as a nexus file
  # open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
  nexus <- readLines(output_name) # open the new nexus file
  ind <- grep("BEGIN CHARACTERS",nexus)+2 # find which line
  nexus[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE;" # replace the line
  writeLines(nexus,output_name) # output the edited nexus file
}

