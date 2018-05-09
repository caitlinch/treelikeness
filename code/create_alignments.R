# Functions to create alignments for simulations with a coalescent or phylogenetic approach

# Create a function to make SimBac alignments
SimBac_make1 <- function(simbac_path, output_folder, ntaxa, nsites, gap, mutation_rate = 0.01, internal_recombination, external_recombination, id = ""){
  # note - site specific mutation rate defaults to 0.01 when not specified
  # Good default gap size is 1,000,000 (1000000)
  output_file <- paste0(output_folder,"SimBac_",ntaxa,"_",nsites,"_",internal_recombination,"_",external_recombination,"_",id,".fasta")
  simbac_command <- paste0(simbac_path," -N ",ntaxa," -B ",nsites," -G ",gap," -T ",mutation_rate," -R ",internal_recombination,
                           " -r ",external_recombination," -o ",output_file)
  system(simbac_command)
}

# Create extreme phylogenetic alignments
# 1. Simulate a tree
# simulate a birth-death tree on a fixed number of extant taxa
# n = number of taxa, numbsim = # of trees to simulate, lambda = speciation rate
# mu = extinction rate, 
tree_sim_depth <- sim.bd.taxa(10, numbsim = 1, lambda = 0.5, mu = 0, frac = 1, complete = FALSE, stochsampling = TRUE)[[1]]
tree_sim_taxa  <- sim.bd.taxa(160, numbsim = 1, lambda = 0.5, mu = 0, frac = 1, complete = FALSE, stochsampling = TRUE)[[1]]
tree_sim_J     <- sim.bd.taxa(10, numbsim = 1, lambda = 0.5, mu = 0, frac = 1, complete = FALSE, stochsampling = TRUE)[[1]]
# 2. Simulate rate variation
# set parameters for creating rate variation in tree
mol_rate <- 0.1
mrate_sd <- 0.1
# make the rate variation in the tree
rate_var_depth <- tree_sim_depth
rate_var_depth$edge.length <- tree_sim_depth$edge.length * rlnorm(length(tree_sim_depth$edge.length), meanlog = log(mol_rate), sdlog = mrate_sd)
rate_var_taxa <- tree_sim_taxa
rate_var_taxa$edge.length <- tree_sim_taxa$edge.length * rlnorm(length(tree_sim_taxa$edge.length), meanlog = log(mol_rate), sdlog = mrate_sd)
rate_var_J <- tree_sim_J
rate_var_J$edge.length <- tree_sim_J$edge.length * rlnorm(length(tree_sim_J$edge.length), meanlog = log(mol_rate), sdlog = mrate_sd)
# scale tree to have total depth of 0.6 million years (CHECK THIS VALUE & FIND A BIOLOGICAL REASON FOR IT)
rate_var_depth <- rescale(rate_var_depth,"depth",5)
rate_var_taxa <- rescale(rate_var_taxa,"depth",0.1)
rate_var_J <- rescale(rate_var_J,"depth",0.1)
# 3. Perform an SPR move for the rate_var_J simulation
rate_var_K <- rSPR(rate_var_J,moves=1)
# 4. Simulate an alignment 
# simulate the sequence alignment based on the tree and number of sites
dna_sim_depth <- as.DNAbin(simSeq(rate_var_depth, l = 1300))
dna_sim_taxa <- as.DNAbin(simSeq(rate_var_taxa, l = 1300))
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

# Create a function to make phylogenetic alignments (as outlined in simulation scheme)
phylo_make1 <- 