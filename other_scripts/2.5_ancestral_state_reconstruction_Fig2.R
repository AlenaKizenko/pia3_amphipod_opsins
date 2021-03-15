#install.packages("phytools")
library(phytools)

## According to the tutorial here: http://www.phytools.org/eqg2015/asr.html 
## or here:
## http://www.phytools.org/Cordoba2017/ex/8/Anc-states-discrete.html

my.tree.text <- readLines("data/2.4_Fig2_tree_rnaspades_nt.nwk")
## https://www.biostars.org/p/294222/ ## how to make the tree rooted
my.tree<-read.tree(text = gsub(";", ":0;", my.tree.text))

## Metadata: the number of opsins. These are discrete values
metadat <- read.csv("data/2.4_metadata.csv")
MWSs <- metadat$MWS
names(MWSs) <- metadat$Species

LWSs <- metadat$LWS
names(LWSs) <- metadat$Species

## example 1 tree
#mtree <- make.simmap(my.tree, MWSs, model="ARD")
#mtree

## but we'll run the analysis multiple times
mtrees <- make.simmap(my.tree, MWSs, nsim = 1000, model="ARD")
pd <- summary(mtrees)

svg("FigS2C.svg")
cols = c("black", "red")
plot(pd, fsize=0.8, cex = 0.4, ftype="i", ylim=c(-2,Ntip(my.tree)))
add.simmap.legend(colors=cols,prompt=FALSE, x=0, y=0, vertical=TRUE, leg = c("zero MWS opsins", "one MWS opsin"), fsize = 0.8)
dev.off()
