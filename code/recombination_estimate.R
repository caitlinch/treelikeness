# R script to estimate the number of recombination events in mammals and in insects
# All values gathered from literature - citations in Draft Simulation Scheme (located on the Dropbox for this project)

# Data gathered from literature
insects_rate <- c(9.7, 19, 16, 4.4, 14.0, 6.2, 4.8, 3.2, 2.5, 5.4, 1.0, 0.3, 0.8, 0.1, 2.6, 0.3, 3.1, 1.6, 
                  6.1, 5.5, 1.2, 2.3, 0.8, 2.9, 4.0, 8.3, 5, 4.1, 6.9) # recombination rate in cM/Mb
insects_rate <- insects_rate*10e6 # insects rate in cM/bp
insects_rate <- insects_rate*10e-2 # insects rate in M/bp <- 1 crossover/generation/bp
insects_geneLength <- c(1346) # mean gene length for eukaryotes in bp
insects_treeDepth <- c(479,434,420,407,396) # origin of insects clade in Mya
insects_genLength <- c(9.37,12.23, 48.12, 51.84, 22.97,55.57, 28.16,30.07) # generation length in days
insects_genLength <- (insects_genLength/(365*10e6)) # generation length in Mya

mammals_rate <- c(0.92,0.6,0.56,1.26,1.04,1.23,1.133,0.734,1.026,1.100,0.762,1.196) # recombination rate in cM/Mb
mammals_rate <- mammals_rate*10e6 # mammals rate in cM/bp
mammals_rate <- mammals_rate*10e-2 # mammals rate in M/bp <- 1 crossover/generation/bp
mammals_geneLength <- c(1346) # mean gene length for eukaryotes in bp
mammals_treeDepth <- c(88,90,80,100,84,122,92.1,108.1,93.9,108.7) # origin of Placentalia clade in Mya
mammals_genLength <- c(880,1700,1095,3410,4290,3438,3190,6200,73,110,730,1095,1095,1460,1460,1825,1825,3650) # generation length in days
mammals_genLength <- (mammals_genLength/(365*10e6)) # generation length in Mya

# Calculate the estimated recombination rate
# To get number of generations: # gens. = (age of clade)/(time for one generation)
# Make sure both generation time and tree depth have units of Mya
insects_numGensUpper <- max(insects_treeDepth)/(min(insects_genLength)) # Get the largest number of generations that could have happened
insects_numGensLower <- min(insects_treeDepth)/(max(insects_genLength)) # Get the smallest number of generations that could have happened
insects_numGensMean   <- mean(insects_treeDepth)/(mean(insects_genLength)) # Get the mean number of generations that could have happened

mammals_numGensUpper <- max(mammals_treeDepth)/(min(mammals_genLength)) # Get the largest number of generations that could have happened
mammals_numGensLower <- min(mammals_treeDepth)/(max(mammals_genLength)) # Get the smallest number of generations that could have happened
mammals_numGensMean   <- mean(mammals_treeDepth)/(mean(mammals_genLength)) # Get the mean number of generations that could have happened

# To get # of crossing over points: # cross. points = crossing over rate / length of gene [((cM/Mb)/Mb) = (cM/Mb)/(Mb/1) = (cM/Mb)*(1/Mb)
insects_crossPointsUpper <- max(insects_rate)/insects_geneLength
insects_crossPointsLower <- min(insects_rate)/insects_geneLength
insects_crossPointsMean <- mean(insects_rate)/insects_geneLength

mammals_crossPointsUpper <- max(mammals_rate)/mammals_geneLength
mammals_crossPointsLower <- min(mammals_rate)/mammals_geneLength
mammals_crossPointsMean <- mean(mammals_rate)/mammals_geneLength

# To get recombination rate: recombination rate = (# cross. points)/(# gens)
insects_minRate <- insects_crossPointsLower/insects_numGensLower
insects_maxRate <- insects_crossPointsUpper/insects_numGensUpper
insects_meanRate <- insects_crossPointsMean/insects_numGensMean

mammals_minRate <- mammals_crossPointsUpper/mammals_numGensUpper
mammals_maxRate <- mammals_crossPointsLower/mammals_numGensLower
mammals_meanRate <- mammals_crossPointsMean/mammals_numGensMean

#To output csv:
type <- c("insects", "insects", "insects", "mammals", "mammals", "mammals")
id <- c("min","mean","max","min","mean","max")
value <- c(insects_minRate,insects_meanRate,insects_maxRate,mammals_minRate,mammals_meanRate,mammals_maxRate)
df <- data.frame(type,id,value)
names(df) <- c("clade","id","recombination_rate")
file <- "/Users/caitlin/Repositories/treelikeness/processed_data/draftSimScheme/recombination_rates.csv"
write.csv(df,file=file)

# To output histograms:
histo <- qplot(insects_rate, geom="histogram",binwidth=2.5e7,main="Histogram for Recombination Rate in Insects",
              xlab="Recombination Frequency (bp/cM)", ylab = "Frequency",fill=I("dark gray"),col=I("black")) + 
  theme(plot.title = element_text(hjust = 0.5))
file_name <- "/Users/caitlin/Repositories/treelikeness/processed_data/draftSimScheme/recombinationRate_insects_histogram.pdf"
ggsave(filename = file_name, plot = histo,device=pdf)

histo <- qplot(mammals_rate, geom="histogram",binwidth = 1e6,main="Histogram for Recombination Rate in Mammals",
                      xlab="Recombination Frequency (bp/cM)", ylab = "Frequency",fill=I("dark gray"),col=I("black")) + 
       theme(plot.title = element_text(hjust = 0.5))
file_name <- "/Users/caitlin/Repositories/treelikeness/processed_data/draftSimScheme/recombinationRate_mammals_histogram.pdf"
ggsave(filename = file_name, plot = histo,device=pdf)

histo <- qplot(insects_genLength*(365*10e6), binwidth=10,geom="histogram",main="Histogram for Generation Length in Insects",
               xlab="Generation Length (days)", ylab = "Frequency",fill=I("dark gray"),col=I("black")) + 
  theme(plot.title = element_text(hjust = 0.5))
file_name <- "/Users/caitlin/Repositories/treelikeness/processed_data/draftSimScheme/generationLength_insects_histogram.pdf"
ggsave(filename = file_name, plot = histo,device=pdf)

histo <- qplot(mammals_genLength*(365*10e6),binwidth = 500, geom="histogram",main="Histogram for Generation Length in Mammals",
                      xlab="Generation Length (days)", ylab = "Frequency",fill=I("dark gray"),col=I("black")) + 
       theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(breaks=seq(0,6500,500))
file_name <- "/Users/caitlin/Repositories/treelikeness/processed_data/draftSimScheme/generationLength_mammals_histogram.pdf"
ggsave(filename = file_name, plot = histo,device=pdf)
