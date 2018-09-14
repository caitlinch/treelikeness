names <- c("08taxa_balanced_LHS.txt", 
           "08taxa_balanced_RHS_nonreciprocal_ancient_1event.txt", 
           "08taxa_balanced_RHS_nonreciprocal_close_1event.txt",
           "08taxa_balanced_RHS_nonreciprocal_divergent_1event.txt", 
           "08taxa_balanced_RHS_reciprocal_ancient_1event.txt", 
           "08taxa_balanced_RHS_reciprocal_close_1event.txt",
           "08taxa_balanced_RHS_reciprocal_divergent_1event.txt",
           "32taxa_balanced_LHS.txt",
           "32_1clade_noevent.txt",
           "32_1clade_nonreciprocalEvent.txt",
           "32_1clade_reciprocalEvent.txt"
           )
out_names_png <- paste0("/Users/caitlincherryh/Documents/Results/plots/seminarPlots/attempt_2/",gsub(".txt",".png",names))
out_names_pdf <- paste0("/Users/caitlincherryh/Documents/Results/plots/seminarPlots/attempt_2/",gsub(".txt",".pdf",names))
names <- paste0("/Users/caitlincherryh/Documents/Repositories/treelikeness/trees/",names)
for (i in 1:length(names)){
  t <- read.tree(names[i])
  png(filename = out_names_png[i])
  plot(t)
  dev.off()
  pdf(file = out_names_pdf[i])
  plot(t)
  dev.off()
}
