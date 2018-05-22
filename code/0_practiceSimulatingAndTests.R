# Code to simulate some alignments, run test statistics on them and see what happens

# Files for SplitsTree practice
alignment_path <- "/Users/caitlincherryh/Documents/test_splitstree/Phylo_20_1300_1_K0.5_tests.nexus"
splitstree_path <- "/Users/caitlincherryh/Documents/Executables/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
iqtree_path       <- "/Applications/iqtree/bin/iqtree" # location of IQ-tree program 

# Open packages
library(TreeSim)
library(phytools)
library(seqinr)
library(ape)
library(phangorn)
library(base)
library(tictoc)
library(ggplot2)
library(reshape2)

# Set working directory
maindir <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/" # for work computer
setwd(maindir)

# Source files for functions
source(paste0(maindir,"code/split_decomposition.R"))
source(paste0(maindir,"code/parametric_bootstrap.R"))
source(paste0(maindir,"code/test_statistic.R"))
source(paste0(maindir,"code/create_alignments.R"))

## If running 3seq remotely, need to specify the ptable before the first run so that 3seq knows where to look
#system("./3seq -f mtDNA.aln -ptable PvalueTable500 -id myFirstRun") # this command associates the Ptable with 3seq - needs a sample alignment to run correctly.

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
for (al in alignments){
  call.IQTREE(iqtree_path,al)
}

# plot trees
test_folder <- "/Users/caitlincherryh/Documents/Repositories/treelikeness/raw_data/testAlignments/"
all_files <- paste0(test_folder,list.files(test_folder))
treefiles <- all_files[grep("treefile",all_files)]
for (treefile in treefiles){
  tree_id <- tail(strsplit(treefile,"/")[[1]],1)
  tree_output_file <- paste0("/Users/caitlincherryh/Documents/TestAlignmentResults/",tree_id,".pdf")
  tree <- read.tree(treefile)
  pdf(file = tree_output_file)
  plot(tree)
  dev.off()
}

setwd(test_folder)
## Run test statistics on these alignments
# Set timer
tic("alignments")
# Set output directory
output_folder <- "/Users/caitlincherryh/Documents/TestAlignmentResults/"
# Set alignment directory
# Open the SimBac alignments
# Set IQ-TREE path
iqtree_path       <- "/Applications/iqtree/bin/iqtree" # location of IQ-tree program 
# Set splitstree path
SplitsTree4_path <- "/Users/caitlincherryh/Documents/Executables/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
# Create storage dataframe
df <- data.frame(matrix(nrow=length(alignments),ncol=8)) # create an empty dataframe of the correct size
string <- c() # empty string to store all the info
row_num <- 1 # make an id for the row number (used to slot the alignment results into the right row in the dataframe)
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
    seq_command <- paste0(seq_path," -f ", al," -d -id ",id)
    system(seq_command) #call 3SEQ
  } else if (filetype == "nexus"){
    # Phipack only reads in Phylip or fasta format - need to convert if the alignment is a nexus file
    data = read.nexus.data(al) # read in nexus format alignment
    fasta.name <- paste0(al,".fasta") # make a name for the fasta alignment by adding .fasta (super original ;) )
    write.fasta(sequences = data,names = names(data), file.out = fasta.name) # output alignment as a fasta format
    phi_command <- paste0(phi_path," -f ",fasta.name) # assemble system command as above
    system(phi_command) # run PHI test on the new fasta alignment
    
    seq_path <- "/Applications/3seq/3seq"
    seq_command <- paste0(seq_path," -f ", fasta.name," -id ",id)
    system(seq_command) #call 3SEQ
  }
  # Extract significance from Phi Pack output
  print("get phi")
  phi_file <- paste0(aldir,"Phi.log")
  phi_file <- readLines(phi_file)
  ind      <- grep("PHI",phi_file)
  phi_sig <- as.numeric(strsplit(phi_file[15],":")[[1]][2])
  
  # Extract output from 3Seq output
  seq_dir <- test_folder # location of 3seq executable
  seq_files <- list.files(seq_dir) # get the list of 3seq output files
  seq_files <- seq_files[grep(id,seq_files)] # prune to only include files for this id
  seq_file <- paste0(seq_dir,seq_files[grep("3s.log",seq_files)]) # get full path to log file
  seq_log <- readLines(seq_file) # open file
  ind      <- grep("Number of recombinant triplets",seq_log) # find the number of recombinant triplets line index
  print("trips")
  num_trips <- seq_log[ind]
  num_trips <- strsplit(num_trips,":")[[1]][2] # extract the number of recombinant triplets
  num_trips <- trimws(num_trips) # trim the whitespace from the number of triplets
  ind      <- grep("Number of distinct recombinant sequences",seq_log) # find the number of distinct recombinant sequences line index
  print("recomb")
  num_dis <- seq_log[ind]
  num_dis <- strsplit(num_dis,":")[[1]][2] # extract the number of distinct recombinant sequences
  num_dis <- trimws(num_dis) # trim the whitespace from the number of distinct recombinant sequences
  # null hypothesis is of clonal evolution - need significant p-value to accept the alternative hypothesis
  ind      <- grep("Rejection of the null hypothesis of clonal evolution",seq_log) # find the p value line index
  print("sig")
  seq_sig <- seq_log[ind]
  seq_sig <- strsplit(seq_sig,"=")[[1]][2] # extract the p value
  seq_sig <- trimws(seq_sig) # trim the whitespace from the number of distinct recombinant sequences
  
  
  # run pdm ratio
  print("pdm")
  pdmr <- pdm.ratio(iqpath = iqtree_path, path = al)
  
  # run normalised.pdm.difference.sum
  print("npds")
  npds <- normalised.pdm.diff.sum(iqpath = iqtree_path, path = al)
  
  # run split decomposition
  print("sd")
  sd <- SplitsTree.decomposition.statistic(iqpath = iqtree_path, splitstree_path = SplitsTree4_path, path = al)
  
  # Collectt results
  row <- c(al,phi_sig,num_trips,num_dis,seq_sig,pdmr,npds,sd) # collect all the information
  print(row)
  df[row_num,] <- row # assign information to correct row in dataframe
  string <- c(string,row) # update string that just contains all data (only doing this in case df doesn't work)
  row_num <- row_num + 1 # iterate up the row number to do the next row in the dataframe
}

# Format output dataframe
names(df) <- c("alignment","PHI","3SEQ_num_recombinant_triplets","3SEQ_num_distinct_recombinant_sequences","3SEQ_p_value","pdm_ratio","pdm_difference","split_decomposition")
write.csv(df,file = "/Users/caitlincherryh/Documents/TestAlignmentResults/test_alignments_results.csv")
toc()

# Make some plots
# Create dataframes for each of the three types of test simulations
phylo_plot_df <- df[1:5,]
K <- c(0,0.01,0.05,0.1,0.5)
phylo_plot_df <- cbind(phylo_plot_df,K)
phylo_plot_df <- phylo_plot_df[,c(9,2,5,6,7,8)]
phylo_plot_df <- melt(phylo_plot_df, id.vars="K")
phylo_plot_df[,2] <- as.character(phylo_plot_df[,2])
phylo_plot <- ggplot(phylo_plot_df,aes(x=K,y=value,col=variable))+geom_point()+
  labs(x = "K (% alignment 2)",  y = "Statistic value", title = "Test Statistic Scores - Phylogenetic")
file_name <- "/Users/caitlincherryh/Documents/TestAlignmentResults/testStatisticScores_phylogenetic.pdf"
ggsave(filename = file_name, plot = phylo_plot, device = "pdf")


internal_plot_df <- df[c(6,10:12),]
internal_recombination <- c(0,0.01,0.1,0.2)
internal_plot_df <- cbind(internal_plot_df,internal_recombination)
internal_plot_df <- internal_plot_df[,c(9,2,5,6,7,8)]
internal_plot_df <- melt(internal_plot_df,id.vars="internal_recombination")
internal_plot_df[,2] <- as.character(internal_plot_df[,2])
internal_plot <- ggplot(internal_plot_df,aes(x=internal_recombination,y=value,col=variable))+geom_point()+
  labs(x = "Internal recombination (%)", y = "Statistic value", title = "Test Statistic Scores - Internal recombination")
file_name <- "/Users/caitlincherryh/Documents/TestAlignmentResults/testStatisticScores_internalRecombination.pdf"
ggsave(filename = file_name, plot = internal_plot, device = "pdf")

external_plot_df <- df[6:9,]
external_recombination <- c(0,0.01,0.05,0.1)
external_plot_df <- cbind(external_plot_df,external_recombination)
external_plot_df <- external_plot_df[,c(9,2,5,6,7,8)]
external_plot_df <- melt(external_plot_df,id.vars="external_recombination")
external_plot_df[,2] <- as.character(external_plot_df[,2])
external_plot <- ggplot(external_plot_df,aes(x=external_recombination,y=value,col=variable))+geom_point()+
  labs(x = "External recombination (%)", y = "Statistic value", title = "Test Statistic Scores - External recombination")
file_name <- "/Users/caitlincherryh/Documents/TestAlignmentResults/testStatisticScores_externalRecombination.pdf"
ggsave(filename = file_name, plot = external_plot, device = "pdf")