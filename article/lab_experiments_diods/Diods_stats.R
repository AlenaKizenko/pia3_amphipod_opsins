library(openxlsx)
library(ggbeeswarm)

## Gmfa for the main figure
diods <- read.delim("diods_summary_Gfa.csv", stringsAsFactors = F)

diods$Experiment <- paste0(diods$Photo, "_", diods$Light, "_", diods$Where)
diods$Type <- paste0(diods$Light, "_", diods$Where, "_", substr(diods$Photo, 7, 10))

diods$Light <- factor(diods$Light, levels = c("K", "B", "G", "Y", "R"))

## for the suppl
ggplot(diods, aes(x=Experiment, y=Xm, col=Light)) + expand_limits(y=c(0, 60)) + 
  geom_boxplot(outlier.shape = NA) + geom_beeswarm() + 
  scale_color_manual(values = c("black", "blue", "green", "orange", "red")) +
  geom_hline(yintercept = c(30), linetype = "dotted") + 
  coord_flip() +
  theme_bw(base_size = 14)

diods[diods$Where=="R", "Xm"] <- 60-diods[diods$Where=="R", "Xm"]
#diods[diods$Light=="K", "Xm"] <- 60-diods[diods$Where=="K", "Xm"]

diods.med <- aggregate(formula=Xm~Experiment+Light, data=diods, FUN="median")

pB <- 
ggplot(diods.med, aes(x=Light, y=Xm, col=Light)) + expand_limits(y=c(0, 60)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = .1) + 
  scale_x_discrete(limits = rev(levels(diods.med$Light)), 
                   labels = c("R"="Red", "Y"="Yellow", "G"="Green", "B"="Blue", "K"="Control (no light)")) + 
  geom_hline(yintercept = c(30), linetype = "dotted") + 
  scale_color_manual(values = c("black", "#00C2F9", "#00D302", "#FFAC38", "#CD022D"), guide=F) +
  coord_flip() + ggtitle("Gm. fasciatus") + #ggtitle("E. cyaneus") +
  theme_bw(base_size = 14) + 
  theme(plot.title = element_text(face = "italic", hjust = .5), axis.text = element_text(size = 12)) +
  ylab("Distance from the edge / diod") + xlab("Light source") + 
  geom_text(data = data.frame(x=1:5, y=rep(58, 5), Light = "K"), aes(x=x, y=y), 
            label = c("*", "**", "*", "*", ""), size = 5, nudge_x = -.1)
pB

ggsave("diods_counts_Gfa.svg", width = 8, height = 3)
ggsave("diods_counts_Gfa.png", width = 8, height = 3)


##pairwise wilcox test?
pairwise.wilcox.test(diods.med$Xm, diods.med$Light)
## with holm: everything significant except for UV (probably not enough replicas)

# data:  diods.med$Xm and diods.med$Light 

#K      B      G      Y     
#B 0.0107 -      -      -     
#  G 0.0023 1.0000 -      -     
#  Y 0.0157 1.0000 1.0000 -     
#  R 0.0107 1.0000 1.0000 1.0000



## Ecy for the main figure
diods <- read.delim("diods_summary_Ecy.tsv", stringsAsFactors = F)

diods$Experiment <- paste0(diods$Photo, "_", diods$Light, "_", diods$Where)
diods$Type <- paste0(diods$Light, "_", diods$Where, "_", substr(diods$Photo, 7, 10))

diods$Light <- factor(diods$Light, levels = c("K", "B", "G", "Y", "R"))

## for the suppl
ggplot(diods, aes(x=Experiment, y=Xm, col=Light)) + expand_limits(y=c(0, 60)) + 
  geom_boxplot(outlier.shape = NA) + geom_beeswarm() + 
  scale_color_manual(values = c("black", "blue", "green", "orange", "red")) +
  geom_hline(yintercept = c(30), linetype = "dotted") + 
  coord_flip() +
  theme_bw(base_size = 14)

diods[diods$Where=="R", "Xm"] <- 60-diods[diods$Where=="R", "Xm"]
#diods[diods$Light=="K", "Xm"] <- 60-diods[diods$Where=="K", "Xm"]

diods.med <- aggregate(formula=Xm~Experiment+Light, data=diods, FUN="median")

##pairwise wilcox test?
pairwise.wilcox.test(diods.med$Xm, diods.med$Light)
## with holm: everything significant except for UV (probably not enough replicas)


pA <- 
  ggplot(diods.med, aes(x=Light, y=Xm, col=Light)) + expand_limits(y=c(0, 60)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = .1) + 
  scale_x_discrete(limits = rev(levels(diods.med$Light)), 
                   labels = c("R"="Red", "Y"="Yellow", "G"="Green", "B"="Blue", "K"="Control (no light)")) + 
  geom_hline(yintercept = c(30), linetype = "dotted") + 
  scale_color_manual(values = c("black", "#00C2F9", "#00D302", "#FFAC38", "#CD022D"), guide=F) +
  coord_flip() + ggtitle("E. cyaneus") + #ggtitle("E. cyaneus") +
  theme_bw(base_size = 14) + 
  theme(plot.title = element_text(face = "italic", hjust = .5), axis.text = element_text(size = 12)) +
  ylab("Distance from the edge / diod") + xlab("Light source") #+ 
#  geom_text(data = data.frame(x=1:5, y=rep(58, 5), Light = "K"), aes(x=x, y=y), 
 #           label = c("*", "**", "*", "*", ""), size = 5, nudge_x = -.1)
pA

ggsave("diods_counts_Ecy.svg", width = 8, height = 3)
ggsave("diods_counts_Ecy.png", width = 8, height = 3)



