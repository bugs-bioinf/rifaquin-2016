source("scripts/plotTree.R")
library(phytools)
library(phangorn)

tru <- read.tree("phylogeny/RAxML_bipartitions.NC_000962.b1")
trr <- midpoint(tru)

plotTree(trr, infoFile="scripts/tipinfo.csv", colourNodesBy ="class", tip.colour.cex=3, lwd=2, tipColours = c("blue","red", "green", "black"), tip.labels=T, tipLabelSize=2.5, offset=0.003, legend=F, infoCols=NA, outputPNG="figures/figure_2.png",w=1000,h=1500,closeDev=F)
edgelabels("Euro American", edge=c(137), adj = c(0.5, -0.5), frame="none",cex=2.5)
edgelabels("East Asian", edge=c(138), adj = c(0.5, -0.5), frame="none",cex=2.5)
edgelabels("East African Indian", edge=c(145), adj = c(0.5, -0.5), frame="none",cex=2.5)
dev.off()

