library(ape)

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
out_names_png <- paste0("/Users/caitlincherryh/Dropbox/Honours/6_paper/figs/trees/nolabels_",gsub(".txt",".png",names))
out_names_pdf <- paste0("/Users/caitlincherryh/Dropbox/Honours/6_paper/figs/trees/nolabels_",gsub(".txt",".pdf",names))
names <- paste0("/Users/caitlincherryh/Documents/Repositories/treelikeness/trees/",names)
for (i in 1:length(names)){
  # plot trees with no labels 
  t <- read.tree(names[i])
  png(filename = out_names_png[i], width = 1000, height = 1000, units = "px")
  plot.phylo(t, edge.width = 8, show.tip.label = FALSE)
  dev.off()
  pdf(file = out_names_pdf[i])
  plot.phylo(t, edge.width = 5, show.tip.label = FALSE)
  dev.off()
}

names <- c("08taxa_balanced_LHS.txt", 
           "08taxa_balanced_RHS_nonreciprocal_ancient_1event.txt", 
           "08taxa_balanced_RHS_nonreciprocal_close_1event.txt",
           "08taxa_balanced_RHS_nonreciprocal_divergent_1event.txt", 
           "08taxa_balanced_RHS_reciprocal_ancient_1event.txt", 
           "08taxa_balanced_RHS_reciprocal_close_1event.txt",
           "08taxa_balanced_RHS_reciprocal_divergent_1event.txt"
)
out_names_png <- paste0("/Users/caitlincherryh/Dropbox/Honours/6_paper/figs/trees/labels_",gsub(".txt",".png",names))
out_names_pdf <- paste0("/Users/caitlincherryh/Dropbox/Honours/6_paper/figs/trees/labels_",gsub(".txt",".pdf",names))
names <- paste0("/Users/caitlincherryh/Documents/Repositories/treelikeness/trees/",names)
for (i in 1:length(names)){
  print(i)
  # plot trees with no labels 
  t <- read.tree(names[i])
  t <- rotateConstr(t,rev(t$tip.label))
  png(filename = out_names_png[i], width = 1000, height = 1000, units = "px")
  plot.phylo(t, edge.width = 8, show.tip.label = TRUE, cex = 6)
  dev.off()
  pdf(file = out_names_pdf[i])
  plot.phylo(t, edge.width = 5, show.tip.label = TRUE, cex = 1.8)
  dev.off()
}

names <- c("32taxa_balanced_LHS.txt")
out_names_png <- paste0("/Users/caitlincherryh/Dropbox/Honours/6_paper/figs/trees/labels_",gsub(".txt",".png",names))
out_names_pdf <- paste0("/Users/caitlincherryh/Dropbox/Honours/6_paper/figs/trees/labels_",gsub(".txt",".pdf",names))
names <- paste0("/Users/caitlincherryh/Documents/Repositories/treelikeness/trees/",names)
for (i in 1:length(names)){
  # plot trees with no labels 
  t <- read.tree(names[i])
  t$tip.label <- paste0("S",32:01)
  png(filename = out_names_png[i], width = 1000, height = 1000, units = "px")
  plot.phylo(t, edge.width = 6, show.tip.label = TRUE, cex = 2.5)
  dev.off()
  pdf(file = out_names_pdf[i])
  plot.phylo(t, edge.width = 3, show.tip.label = TRUE, cex = 1.3)
  dev.off()
}

names <- c("32_1clade_noevent.txt","32_1clade_nonreciprocalEvent.txt","32_1clade_reciprocalEvent.txt")
out_names_png <- paste0("/Users/caitlincherryh/Dropbox/Honours/6_paper/figs/trees/labels_",gsub(".txt",".png",names))
out_names_pdf <- paste0("/Users/caitlincherryh/Dropbox/Honours/6_paper/figs/trees/labels_",gsub(".txt",".pdf",names))
names <- paste0("/Users/caitlincherryh/Documents/Repositories/treelikeness/trees/",names)
for (i in 1:length(names)){
  # plot trees with no labels 
  t <- read.tree(names[i])
  png(filename = out_names_png[i], width = 1000, height = 1000, units = "px")
  plot.phylo(t, edge.width = 6, show.tip.label = TRUE, cex = 5)
  dev.off()
  pdf(file = out_names_pdf[i])
  plot.phylo(t, edge.width = 3, show.tip.label = TRUE, cex = 3)
  dev.off()
}
