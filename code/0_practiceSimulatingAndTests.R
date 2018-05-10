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
source(paste0(maindir,"code/create_alignments.R"))

# Create alignments
simbac_path <- "/Users/caitlincherryh/Documents/Repositories/SimBac/SimBac"
output_folder <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments/"
ntaxa <- 20
nsites <- 1300
gap <- 1000000
mutation_rate <- 0.01
id <- "tests"
SimBac.make1(simbac_path, output_folder, ntaxa, nsites, gap, mutation_rate, 0.00, 0.00, id) # Order: mutation_rate, internal recombination, external recombination
SimBac.make1(simbac_path, output_folder, ntaxa, nsites, gap, mutation_rate, 0.00, 0.01, id) # Order: mutation_rate, internal recombination, external recombination
SimBac.make1(simbac_path, output_folder, ntaxa, nsites, gap, mutation_rate, 0.00, 0.05, id) # Order: mutation_rate, internal recombination, external recombination
SimBac.make1(simbac_path, output_folder, ntaxa, nsites, gap, mutation_rate, 0.00, 0.10, id) # Order: mutation_rate, internal recombination, external recombination
SimBac.make1(simbac_path, output_folder, ntaxa, nsites, gap, mutation_rate, 0.01, 0.00, id) # Order: mutation_rate, internal recombination, external recombination
SimBac.make1(simbac_path, output_folder, ntaxa, nsites, gap, mutation_rate, 0.10, 0.00, id) # Order: mutation_rate, internal recombination, external recombination
SimBac.make1(simbac_path, output_folder, ntaxa, nsites, gap, mutation_rate, 0.20, 0.00, id) # Order: mutation_rate, internal recombination, external recombination

# Add extra parameters needed to create the phylogenetic alignments
birth_rate = 0.5
death_rate = 0
tree_age = 1
mol_rate = 0.1
mol_rate_sd = 0.1
K_vector = c(0,0.01,0.05,0.1,0.5) # proportion of second tree in the mosaic alignment (here 0 - 50% in 1% increments for a total of 51 alignments. Must be in decimals)
phylo.make1(output_folder, ntaxa, nsites, birth_rate, death_rate, tree_age, mol_rate, mol_rate_sd, K_vector,id)

# Get a list of all the alignment files
aldir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments/"
setwd(aldir)
alignments <- list.files(aldir) 

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
for (treefile in treefiles){
  tree_id <- strsplit(tail(strsplit(treefile,"/")[[1]],1),"\\.")[[1]][1]
  tree_output_file <- paste0("/Users/caitlincherryh/Documents/TestAlignmentResults/",tree_id,".pdf")
  tree <- read.tree(treefile)
  pdf(tree_output_file)
  plot(tree)
  dev.off()
}

## Run test statistics on these alignments
# Set timer
tic("alignments")
# Set output directory
output_folder <- "/Users/caitlincherryh/Documents/TestAlignmentResults/"
# Set alignment directory
# Open the SimBac alignments
# Set IQ-TREE path
iqtree_path       <- "/Applications/iqtree/bin/iqtree" # location of IQ-tree program 
# Create storage dataframe
df <- data.frame(matrix(nrow=0,ncol=8))
# Run test statistics on each alignment
# Record values for test statistics
for (al in alignments){
  # run PHIPACK and 3seq
  phi_path <- "/Applications/PhiPack/Phi"
  filetype = tail(strsplit(al,"\\.")[[1]],n=1) # extract file format
  id <- tail(strsplit(al,"/")[[1]],n=1)
  if (filetype == "fasta"){
    # if the alignment is already in fasta format, run PhiPack through R
    phi_command <- paste0(phi_path," -f ",al) # assemble system command
    system(phi_command) #call phipack
    
    seq_path <- "/Applications/3seq/3seq"
    seq_command <- paste0(seq_path," -f ", al," -d -p -id ",id)
    system(seq_command) #call 3SEQ
  } else if (filetype == "nexus"){
    # Phipack only reads in Phylip or fasta format - need to convert if the alignment is a nexus file
    data = read.nexus.data(al) # read in nexus format alignment
    fasta.name <- paste0(al,".fasta") # make a name for the fasta alignment by adding .fasta (super original ;) )
    write.fasta(sequences = data,names = names(data), file.out = fasta.name) # output alignment as a fasta format
    phi_command <- paste0(phi_path," -f ",fasta.name) # assemble system command as above
    system(phi_command) # run PHI test on the new fasta alignment
    
    seq_path <- "/Applications/3seq/3seq"
    seq_command <- paste0(seq_path," -f ", fasta.name," -p -id ",id)
    system(seq_command) #call 3SEQ
  }
  # Extract significance from Phi Pack output
  phi_file <- paste0(aldir,"Phi.log")
  phi_file <- readLines(phi_file)
  ind      <- grep("PHI",phi_file)
  phi_sig <- as.numeric(strsplit(phi_file[15],":")[[1]][2])
  
  # Extract output from 3Seq output
  seq_dir <- "/Applications/3seq/" # location of 3seq executable
  seq_files <- list.files(seq_dir) # get the list of 3seq output files
  seq_files <- seq_files[grep(id,seq_files)] # prune to only include files for this id
  seq_file <- paste0(seq_dir,seq_files[grep("log",seq_files)]) # get full path to log file
  seq_log <- readLines(seq_file) # open file
  ind      <- grep("Number of recombinant triplets",seq_log) # find the number of recombinant triplets line index
  num_trips <- seq_log[ind]
  num_trips <- strsplit(num_trips,":")[[1]][2] # extract the number of recombinant triplets
  num_trips <- trimws(num_trips) # trim the whitespace from the number of triplets
  ind      <- grep("Number of distinct recombinant sequences",seq_log) # find the number of distinct recombinant sequences line index
  num_dis <- seq_log[ind]
  num_dis <- strsplit(num_dis,":")[[1]][2] # extract the number of distinct recombinant sequences
  num_dis <- trimws(num_dis) # trim the whitespace from the number of distinct recombinant sequences
  # null hypothesis is of clonal evolution - need significant p-value to accept the alternative hypothesis
  ind      <- grep("Rejection of the null hypothesis of clonal evolution",seq_log) # find the p value line index
  seq_sig <- seq_log[ind]
  seq_sig <- strsplit(seq_sig,"=")[[1]][2] # extract the p value
  seq_sig <- trimws(seq_sig) # trim the whitespace from the number of distinct recombinant sequences
  
  
  # run pdm ratio
  pdmr <- pdm.ratio(iqpath = iqtree_path, path = al)
  
  # run normalised.pdm.difference.sum
  npds <- normalised.pdm.diff.sum(iqpath = iqtree_path, path = al)
  
  # run split decomposition
  sd <- split.decomposition.statistic(iq_path = iqtree_path, path = al)
  
  # Collectt results
  row <- c(al,phi_sig,num_trips,num_dis,seq_sig,pdmr,npds,sd)
  df <- rbind(df,row)
}

# Format output dataframe
names(df) <- c("alignment","PHI","3SEQ_num_recombinant_triplets","3SEQ_num_distinct_recombinant_sequences","3SEQ_p_value","pdm_ratio","pdm_difference","split_decomposition")
write.csv(df,file = "/Users/caitlincherryh/Documents/TestAlignmentResults/test_alignments_results.csv")
toc()

# Make some plots