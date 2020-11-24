## branch length analysis

library(ggtree)

#read the tree
tr <- read.iqtree("data/2.2_allopsins.trim.aln.treefile")

#only edges that contain final leaves (1:149)
edges.keep <- tr@phylo$edge[,2] < 150

hist(tr@phylo$edge.length[edges.keep], breaks = 150)
## B. taurus is around 3.89
abline(v = 0.5, lty=3)

## Finally, the tree (Fig. S2A) was visualized with iTOL (https://itol.embl.de/).
