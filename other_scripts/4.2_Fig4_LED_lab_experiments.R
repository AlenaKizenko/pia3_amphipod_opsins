library(openxlsx)
library(ggbeeswarm)
library(ggpubr)

## Gmfa for the main figure
## read the data
diodsG <- read.delim("data/4.1_summary_Gfa_LED.csv", stringsAsFactors = F)

## Add new variables for the groups (useful for plotting)
diodsG$Experiment <- paste0(diodsG$Photo, "_", diodsG$Light, "_", diodsG$Where)
diodsG$Type <- paste0(diodsG$Light, "_", diodsG$Where, "_", substr(diodsG$Photo, 7, 10))
diodsG$Light <- factor(diodsG$Light, levels = c("K", "B", "G", "Y", "R"))

## A way to plot each experiment as one panel
#ggplot(diodsG, aes(x=Experiment, y=Xm, col=Light)) + expand_limits(y=c(0, 60)) + 
#  geom_boxplot(outlier.shape = NA) + geom_beeswarm() + 
#  scale_color_manual(values = c("black", "blue", "green", "orange", "red")) +
#  geom_hline(yintercept = c(30), linetype = "dotted") + 
#  coord_flip() +
#  theme_bw(base_size = 14)

## Correct for the cases where the light sources was on the right wall of the tank
diodsG[diodsG$Where=="R", "Xm"] <- 60-diodsG[diodsG$Where=="R", "Xm"]

## Compute medians to compare and plot
diodsG.med <- aggregate(formula=Xm~Experiment+Light+Lab.acclimation+Where, data=diodsG, FUN="median")

## Test the significance of the differences
pairwise.wilcox.test(diodsG.med$Xm, diodsG.med$Light)
# data:  diods.med$Xm and diods.med$Light 
#K      B      G      Y     
#B 0.0107 -      -      -     
#  G 0.0023 1.0000 -      -     
#  Y 0.0157 1.0000 1.0000 -     
#  R 0.0107 1.0000 1.0000 1.0000


## And plot
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

## Ecy for the main figure
## The same code as before; added purely for reproducibility
diodsE <- read.delim("data/4.1_summary_Ecy_LED.csv", stringsAsFactors = F)

## Add new variables for the groups (useful for plotting)
diodsE$Experiment <- paste0(diodsE$Photo, "_", diodsE$Light, "_", diodsE$Where)
diodsE$Type <- paste0(diodsE$Light, "_", diodsE$Where, "_", substr(diodsE$Photo, 7, 10))
diodsE$Light <- factor(diodsE$Light, levels = c("K", "B", "G", "Y", "R"))

## A way to plot each experiment as one panel
#ggplot(diodsE, aes(x=Experiment, y=Xm, col=Light)) + expand_limits(y=c(0, 60)) + 
#  geom_boxplot(outlier.shape = NA) + geom_beeswarm() + 
#  scale_color_manual(values = c("black", "blue", "green", "orange", "red")) +
#  geom_hline(yintercept = c(30), linetype = "dotted") + 
#  coord_flip() +
#  theme_bw(base_size = 14)

## Correct for the cases where the light sources was on the right wall of the tank
diodsE[diodsE$Where=="R", "Xm"] <- 60-diodsE[diodsE$Where=="R", "Xm"]

## Compute medians to compare and plot
diodsE.med <- aggregate(formula=Xm~Experiment+Light+Lab.acclimation+Where, data=diodsE, FUN="median")

## Test the significance of the differences
pairwise.wilcox.test(diodsE.med$Xm, diodsE.med$Light)
#     K       B       G       Y      
#  B 0.00016 -       -       -      
#  G 0.00016 1.00000 -       -      
#  Y 0.00016 0.03497 0.24895 -      
#  R 0.00016 0.00093 0.11251 1.00000

## And plot
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


## Finally save the figure
ggarrange(pA, pB + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + xlab(""), ncol = 2, widths = c(1, 0.7), labels = c("A", "B"))
ggsave("Ecy_Gmfa_diods.png", width = 8, height = 3)

### A few more things to know
## Do the acclimation and position of the diod (left/right) matter?
summary(glm(data=diodsG.med, Xm~Light+Lab.acclimation+Where))
summary(glm(data=diodsE.med, Xm~Light+Lab.acclimation+Where))
## No, they don't.