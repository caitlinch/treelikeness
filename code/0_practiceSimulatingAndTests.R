# Code to simulate some alignments, run test statistics on them and see what happens

# Open packages
library(TreeSim)
library(phytools)
library(seqinr)
library(ape)
library(phangorn)
library(base)
library(tictoc)

# Set working directory
maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # for work computer
setwd(maindir)

# Source files for functions
source(paste0(maindir,"code/split_decomposition.R"))
source(paste0(maindir,"code/parametric_bootstrap.R"))
source(paste0(maindir,"code/test_statistic.R"))

# Create alignments
# Run SimBac to get extreme alignments
simbac_path <- "/Users/caitlincherryh/Documents/Repositories/SimBac/Simbac"
test_folder <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments/"
# ./SimBac -N 100 -B 1300 -G 1000000 -R 0.00 -o genomes.fasta -c clonal_frame.nwk
system(paste0(simbac_path," -N 160 -B 1300 -G 1000000 -R 0.00 -o ", test_folder,"simbac_160taxa.fasta -c ",test_folder,"clonal_frame_160taxa.nwk"))
system(paste0(simbac_path," -N 10 -B 1300 -G 1000000 -R 0.2 -o ", test_folder,"simbac_0.2R_int.fasta -c ",test_folder,"clonal_frame_0.2R_int.nwk"))
system(paste0(simbac_path," -N 10 -B 1300 -G 1000000 -r 0.02 -R 0.00 -o ",test_folder, "simbac_0.2R_ext.fasta -c ",test_folder,"clonal_frame_0.2R_ext.nwk"))
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

# Run IQ-tree in each alignment
# Get the list of files in the folder
test_folder <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments/"
alignments <- paste0(test_folder,list.files(test_folder))
# Set IQ-TREE path
iqtree_path       <- "/Applications/iqtree/bin/iqtree" # location of IQ-tree program 
for (al in alignments[1:3]){
  call.IQTREE(iqtree_path,al)
}

# plot trees
test_folder <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments/"
alignments <- paste0(test_folder,list.files(test_folder))
treefiles <- alignments[grep("treefile",alignments)]

## Run test statistics on these alignments
# Set timer
tic("alignments")
# Set output directory
output_folder <- "/Users/caitlincherryh/Documents/TestAlignmentResults/"
# Set alignment directory
aldir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments/"
setwd(aldir)
als <- c("alignment_phylo_5050.nexus","alignment_phylo_depth.nexus","alignment_phylo_taxa.nexus","simbac_0.2R_ext.fasta","simbac_0.2R_int.fasta","simbac_160taxa.fasta")
# Open the SimBac alignments
alignments <- paste0(aldir,als)
# Set IQ-TREE path
iqtree_path       <- "/Applications/iqtree/bin/iqtree" # location of IQ-tree program 
# Create storage dataframe
df <- data.frame(matrix(nrow=0,ncol=6))
# Run test statistics on each alignment
# Record values for test statistics
for (al in alignments){
  # run PHIPACK
  phi_path <- "/Applications/PhiPack/Phi"
  filetype = tail(strsplit(al,"\\.")[[1]],n=1) # extract file format
  if (filetype == "fasta"){
    # if the alignment is already in fasta format, run PhiPack through R
    phi_command <- paste0(phi_path," -f ",al) # assemble system command
    system(phi_command) #call phipack
  } else if (filetype == "nexus"){
    # Phipack only reads in Phylip or fasta format - need to convert if the alignment is a nexus file
    data = read.nexus.data(al) # read in nexus format alignment
    fasta.name <- paste0(al,".fasta") # make a name for the fasta alignment by adding .fasta (super original ;) )
    write.fasta(sequences = data,names = names(data), file.out = fasta.name) # output alignment as a fasta format
    phi_command <- paste0(phi_path," -f ",fasta.name) # assemble system command as above
    system(phi_command) # run PHI test on the new fasta alignment
  }
  phi_file <- paste0(aldir,"Phi.log")
  phi_file <- readLines(phi_file)
  ind      <- grep("PHI",phi_file)
  phi_sig <- as.numeric(strsplit(phi_file[15],":")[[1]][2])
  
  # for now, don't run 3SEQ
  #seq_path <- "/Applications/3seq/3seq"
  #seq_command <- paste0(seq_path," -f ", al," -d -p")
  #system(seq_command) #call 3SEQ
  seq_sig <- 0 # put results from 3SEQ here
  
  # run pdm ratio
  pdmr <- pdm.ratio(iqpath = iqtree_path, path = al)
  
  # run normalised.pdm.difference.sum
  npds <- normalised.pdm.diff.sum(iqpath = iqtree_path, path = al)
  
  # run split decomposition
  sd <- split.decomposition.statistic(iq_path = iqtree_path, path = al)
  
  # Collectt results
  row <- c(al,phi_sig,seq_sig,pdmr,npds,sd)
  df <- rbind(df,row)
}

# Format output dataframe
names(df) <- c("alignment","PHI","3SEQ","pdm_ratio","pdm_difference","split_decomposition")
toc()

# Make some plots