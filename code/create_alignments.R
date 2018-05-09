# Functions to create alignments for simulations with a coalescent or phylogenetic approach

# Create a function to make SimBac alignments
SimBac_make1 <- function(simbac_path, output_folder, ntaxa, nsites, gap, mutation_rate = 0.01, internal_recombination, external_recombination, id = ""){
  # note - site specific mutation rate defaults to 0.01 when not specified in SimBac
  # Good default gap size is 1,000,000 (1000000)
  output_file <- paste0(output_folder,"SimBac_",ntaxa,"_",nsites,"_",internal_recombination,"_",external_recombination,"_",id,".fasta")
  simbac_command <- paste0(simbac_path," -N ",ntaxa," -B ",nsites," -G ",gap," -T ",mutation_rate," -R ",internal_recombination,
                           " -r ",external_recombination," -o ",output_file)
  system(simbac_command)
}

# want to feed in a vector of J % that you want to test

# Create a function to make phylogenetic alignments (as outlined in simulation scheme)
# Provide J_vector in decimals (e.g. 1% = 0.01, 50% = 0.5)
# J is the proportion of the SECOND tree that will be included
phylo_make1 <- function(ntaxa, nsites, birth_rate = 0.5, death_rate = 0, tree_age, mol_rate = 0.1, mol_rate_sd = 0.1, J_vector = c()){
  # 1. Simulate a tree
  # simulate a birth-death tree on a fixed number of extant taxa
  # n = number of taxa, numbsim = # of trees to simulate, lambda = speciation rate [good default = 0.5 from Duchenne and Lanfear (2015)]
  # mu = extinction rate (default for this project = 0)
  tree_sim <- sim.bd.taxa.age(n = ntaxa, numbsim = 1, lambda = birth_rate, mu = death_rate, frac = 1, age = tree_age, mrca = FALSE)[[1]]
  tree_sim$edge.length <- tree_sim$edge.length * (tree_age / max(branching.times(tree_sim)))
  # 2. Simulate rate variation
  # Default for mol_rate and mol_rate_sd = 0.1 as in Duchenne and Lanfear (2015)
  phylo_sim <- tree_sim
  phylo_sim$edge.length <- tree_sim$edge.length * rlnorm(length(tree_sim$edge.length), meanlog = log(mol_rate), sdlog = mol_rate_sd)
  # scale tree to have a total depth of tree age
  phylo_sim <- rescale(phylo_sim,"depth",tree_age)
  # If the J vector is empty, simply simulate DNA along the tree
  if (length(J_vec)==0){
    dna_sim <- as.DNAbin(simSeq(phylo_sim),l = nsites)
  }
  else {
    # If there are elements in the J vector, need to create a 2nd alignment to concatenate at those intervals
    phylo_sim_2 <- rSPR(phylo_sim, moves=1) # perform a single SPR move at random
    K_vector <- 1 - J_vector # proportion of first tree that will be included
    dna_sim_1 <- as.DNAbin(simSeq(phylo_sim),l = K_vector*nsites) # simulate along the entire first tree
    dna_sim_2 <- as.DNAbin(simSeq(phylo_sim_2),l = J_vector*nsites) # simulate along the entire second tree
    
  }
}


# 4. Simulate an alignment 
# simulate the sequence alignment based on the tree and number of sites
dna_sim_depth <- as.DNAbin(simSeq(rate_var_depth, l = 1300))

dna_sim_J <- as.DNAbin(simSeq(rate_var_J, l = 650))
dna_sim_K <- as.DNAbin(simSeq(rate_var_K, l = 650))
# concatenate alignments J and K
dna_sim_5050 <- cbind(dna_sim_J,dna_sim_K)
# save alignment as a nexus file
op_file <- paste0(test_folder,"/alignment_phylo_depth.nexus")
write.nexus.data(dna_sim_depth, file = op_file)
op_file <- paste0(test_folder,"/alignment_phylo_taxa.nexus")
write.nexus.data(dna_sim_taxa, file = op_file)
op_file <- paste0(test_folder,"/alignment_phylo_5050.nexus")
write.nexus.data(dna_sim_5050, file = op_file)
