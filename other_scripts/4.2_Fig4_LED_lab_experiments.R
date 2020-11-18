library(openxlsx)
library(ggbeeswarm)

## Gmfa for the main figure
diodsG <- read.delim("summary_Gfa_LED.tsv.csv", stringsAsFactors = F)

diodsG$Experiment <- paste0(diodsG$Photo, "_", diodsG$Light, "_", diodsG$Where)
diodsG$Type <- paste0(diodsG$Light, "_", diodsG$Where, "_", substr(diodsG$Photo, 7, 10))

diodsG$Light <- factor(diodsG$Light, levels = c("K", "B", "G", "Y", "R"))

## for the suppl
ggplot(diodsG, aes(x=Experiment, y=Xm, col=Light)) + expand_limits(y=c(0, 60)) + 
  geom_boxplot(outlier.shape = NA) + geom_beeswarm() + 
  scale_color_manual(values = c("black", "blue", "green", "orange", "red")) +
  geom_hline(yintercept = c(30), linetype = "dotted") + 
  coord_flip() +
  theme_bw(base_size = 14)

diodsG[diodsG$Where=="R", "Xm"] <- 60-diodsG[diodsG$Where=="R", "Xm"]
#diodsG[diodsG$Light=="K", "Xm"] <- 60-diodsG[diodsG$Where=="K", "Xm"]

#diodsG.med <- aggregate(formula=Xm~Experiment+Light, data=diodsG, FUN="median")
diodsG.med <- aggregate(formula=Xm~Experiment+Light+Lab.acclimation+Where, data=diodsG, FUN="median")


##pairwise wilcox test?
pairwise.wilcox.test(diodsG.med$Xm, diodsG.med$Light)
## with holm: everything significant except for UV (probably not enough replicas)

# data:  diods.med$Xm and diods.med$Light 

#K      B      G      Y     
#B 0.0107 -      -      -     
#  G 0.0023 1.0000 -      -     
#  Y 0.0157 1.0000 1.0000 -     
#  R 0.0107 1.0000 1.0000 1.0000



pB <- 
ggplot(diodsG.med, aes(x=Light, y=Xm, col=Light)) + expand_limits(y=c(0, 60)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = .1) +#, aes(shape=Lab.acclimation)) +  or shape=Where
  scale_x_discrete(limits = rev(levels(diodsG.med$Light)), 
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

# ggsave("diods_counts_Gfa.svg", width = 8, height = 3)
# ggsave("diods_counts_Gfa.png", width = 8, height = 3)



## Ecy for the main figure
diodsE <- read.delim("summary_Ecy_LED.csv", stringsAsFactors = F)

diodsE$Experiment <- paste0(diodsE$Photo, "_", diodsE$Light, "_", diodsE$Where)
diodsE$Type <- paste0(diodsE$Light, "_", diodsE$Where, "_", substr(diodsE$Photo, 7, 10))

diodsE$Light <- factor(diodsE$Light, levels = c("K", "B", "G", "Y", "R"))

## for the suppl
ggplot(diodsE, aes(x=Experiment, y=Xm, col=Light)) + expand_limits(y=c(0, 60)) + 
  geom_boxplot(outlier.shape = NA) + geom_beeswarm() + 
  scale_color_manual(values = c("black", "blue", "green", "orange", "red")) +
  geom_hline(yintercept = c(30), linetype = "dotted") + 
  coord_flip() +
  theme_bw(base_size = 14)

diodsE[diodsE$Where=="R", "Xm"] <- 60-diodsE[diodsE$Where=="R", "Xm"]
#diodsE[diodsE$Light=="K", "Xm"] <- 60-diodsE[diodsE$Where=="K", "Xm"]

diodsE.med <- aggregate(formula=Xm~Experiment+Light+Lab.acclimation+Where, data=diodsE, FUN="median")

##pairwise wilcox test?
pairwise.wilcox.test(diodsE.med$Xm, diodsE.med$Light)
## with holm: everything significant except for UV (probably not enough replicas)


pA <- 
  ggplot(diodsE.med, aes(x=Light, y=Xm, col=Light)) + expand_limits(y=c(0, 60)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = .1) +#, aes(shape=Lab.acclimation)) +  or shape = Where
  scale_x_discrete(limits = rev(levels(diodsE.med$Light)), 
                   labels = c("R"="Red", "Y"="Yellow", "G"="Green", "B"="Blue", "K"="Control (no light)")) + 
  geom_hline(yintercept = c(30), linetype = "dotted") + 
  scale_color_manual(values = c("black", "#00C2F9", "#00D302", "#FFAC38", "#CD022D"), guide=F) +
  coord_flip() + ggtitle("E. cyaneus") + #ggtitle("E. cyaneus") +
  theme_bw(base_size = 14) + 
  theme(plot.title = element_text(face = "italic", hjust = .5), axis.text = element_text(size = 12)) +
  ylab("Distance from the edge / diod") + xlab("Light source") + 
  geom_text(data = data.frame(x=1:5, y=rep(58, 5), Light = "K"), aes(x=x, y=y), 
           label = c("***", "***", "***", "***", ""), size = 5, nudge_x = -.1)
pA

# ggsave("diods_counts_Ecy.svg", width = 8, height = 3)
# ggsave("diods_counts_Ecy.png", width = 8, height = 3)


##pairwise wilcox test?
pairwise.wilcox.test(diods.med$Xm, diods.med$Light)

# data:  diods.med$Xm and diods.med$Light 

# K       B       G       Y      
# B 9e-06   -       -       -      
#   G 9e-06   1.00000 -       -      
#   Y 9e-06   0.03497 0.24895 -      
#   R 0.00019 0.00093 0.11251 1.00000
# 
# P value adjustment method: holm 


### does acclimation matter?
summary(glm(data=diodsG.med, Xm~Light+Lab.acclimation+Where))

summary(glm(data=diodsE.med, Xm~Light+Lab.acclimation+Where))


library(ggpubr)
#ggarrange(pA, pB, nrow = 2, labels = c("A", "B"))
ggarrange(pA, pB + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + xlab(""), ncol = 2, widths = c(1, 0.7), labels = c("A", "B"))
ggsave("Ecy_Gmfa_diods.png", width = 8, height = 3)
