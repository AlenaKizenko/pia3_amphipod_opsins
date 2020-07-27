options(stringsAsFactors = F)
library(ggtree)
library(ggplot2)
## if installing ggtree doesn't work
#install.packages("tidytree") #it is on CRAN!

## install the main package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree", version = "3.8")

#tr <- read.iqtree("tree_rnaspades_nt.nwk")
tr <- read.iqtree("tree_rnaspades_aa.nwk")

label <- tr@phylo$node.label
alrt <- as.numeric(sub("/.*", "", label))
bigalrt <- alrt > 70 & !is.na(alrt)
##the first number
bayes <- as.numeric(sub(".*/", "", label))
bigbayes <- bayes > 0.7 & !is.na(bayes)
## subset
newlabel <- ifelse(bigalrt & bigbayes, "*", "")
#newlabel <- ifelse(bigalrt & bigbayes, label, "")
tr@phylo$node.label <- newlabel
#tr@phylo$tip.label <- gsub("_", " ", tr@phylo$tip.label, )


ptr <- ggtree(tr) + geom_tiplab()


meta.df <- read.csv("metadata.csv", row.names = 2)[, 2:3]
#meta.df[,2] <- factor(meta.df[,2])
#meta.df[,1] <- factor(meta.df[,1])

ptr2 <- gheatmap(ptr, meta.df, width=.2, offset=.2, colnames=T, color = "black", low = "white", high = "black") %>% 
  scale_x_ggtree

ptr2 + geom_nodelab(size = 7, color = "darkred", nudge_x = -0.005, nudge_y = -0.3, geom="text")


#ggsave("nt_tree.svg")
ggsave("aa_tree.svg")

#+  geom_text2(data = ptr$data[(!ptr$data$isTip), ], 
#                   aes(subset=!isTip, label=tr@phylo$node.label), size = 1, 
#                   family = "arial", col = "darkred") 
#+ geom_nodelab(size = 3, nudge_x = -0.0001, nudge_y = 0.0, geom="text")
#+ geom_nodepoint(aes(subset = as.logical(newlabel == "")))


# ptr2 + scale_fill_brewer(palette="Set2") + theme_tree2() + 
#   scale_y_continuous(expand=c(0, 0.6)) + xlab("Time") +
#   theme(legend.text=element_text(size=8), 
#          legend.key.height=unit(.5, "cm"),
#          legend.key.width=unit(.4, "cm"), 
#          legend.position=c(.13, y=.945),
#          axis.text.x=element_text(size=10), 
#          axis.title.x = element_text(size=12))

library(RColorBrewer)


p1 <- ggtree(tr) + 
#  geom_hilight(441, fill = "#24ff24", extendto=0.7) +
#  geom_hilight(656, fill = "#fd6cb5", extendto=0.64) +
#  geom_hilight(679, fill = "#6ab5fa", extendto=0.6) +
#  geom_hilight(833, fill = "#ffeb8a", extendto=0.89) +
  geom_nodelab(size = 2, nudge_x = -0.003, nudge_y = 0.2, geom="text") +
  #geom_nodepoint(aes(subset = as.logical(newlabel == ""))) + 
  #geom_point2(aes(subset = as.logical(newlabel == "") & !isTip)) + 
  #geom_nodepoint(aes(subset = 
  #          (as.numeric(sub("/.*", "", label)) > 70 | 
  #            as.numeric(sub(".*/", "", label) > 0.7))  & !isTip)) + 
#  geom_tiplab(size = 2.5, fontface = "italic", align = T) +
  

p1
#identify(p1)

#ggsave("tree.svg")


p1 <- ggtree(tr) + 
  geom_nodelab(size = 2, nudge_x = -0.003, nudge_y = 0.2, geom="text") +
  geom_tiplab(size = 3.5, fontface = "italic") + 
  ggplot2::xlim(0, 1)

p1  

## with phyloseq
## good idea... but not quite what we wanted
### https://yulab-smu.github.io/treedata-book/chapter9.html#ggtree-for-other-tree-like-structure

# numopsins <- t(meta.df[,2:3])
# colnames(numopsins) <- meta.df[,1]
# 
# meta.df2 <- meta.df[,2:3] 
# row.names(meta.df2) <- meta.df[,1]
# 
# otut <- otu_table(meta.df2, taxa_are_rows = T)
# trr <- read.tree("tree_rnaspades_nt.nwk")
# samd <- sample_data(data.frame(SampleID = sample_names(otut), colID = c("blue", "green"), row.names = sample_names(otut)))
# 
# ps <- phyloseq(otut, trr, samd)
# 
# ggtree(ps) +
#   geom_nodelab(size = 2, nudge_x = -0.003, nudge_y = 0.2, geom="text") +
#   geom_tiplab(align = T, size = 3.5, fontface = "italic") + 
#   ggplot2::xlim(0, 1) +
#   geom_point2(aes(subset=isTip, x=x+hjust, color = colID), alpha=.3, na.rm=T) + 
#   theme(legend.position="right")

