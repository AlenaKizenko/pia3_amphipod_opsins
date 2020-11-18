library(reshape2)
library(scales)
library(ggplot2)
library(ggbeeswarm)

## data input (it's Table S2)
aln.data <- read.csv("alignment_to_Gpu_Gammaridae.csv")

aln.data <- droplevels(aln.data[aln.data$Include.to.main.fig,])

# sum(aln.data[aln.data$Species == "Eulimnogammarus cyaneus", "M.Seqs"])
# sum(aln.data[aln.data$Species == "Eulimnogammarus verrucosus", "M.Seqs"])
# sum(aln.data[aln.data$Species == "Gammarus lacustris", "M.Seqs"])
# 
# sum(aln.data[aln.data$Species == "Gammarus pulex", "M.Seqs"])


# species <- read.csv("../species_tree/metadata.csv")[, "Species"]
# species <- gsub("_", " ", species)
# species <- c(species, "Eulimnogammarus verrucosus")


aln.data.show <- melt(aln.data[, c("Group", "M.Seqs", "MWS.RPM", "LWS1.RPM", "LWS2.RPM")],
                      variable.name = "opstype", value.name = "RPM", id.vars = c("Group", "M.Seqs"))



# aln.data.show <- melt(aln.data[aln.data$Species %in% species, c("Species", "MWS.RPM", "LWS.RPM")],
#                        variable.name = "opstype", value.name = "RPM")

#aln.data.show <- aln.data.show[aln.data.show$Species %in% ]

# ggplot(aln.data.show) +
#   geom_boxplot(aes(y = RPM, x = Species), varwidth = T) +
#   geom_point(aes(y = RPM, x = Species), size = .1) +
#   #geom_beeswarm(aes(y = RPM, x = Species)) + 
#   coord_flip() + 
#   facet_wrap( ~ opstype, scales = "free_x") +
#   theme_light(base_size = 24)

#aggregate(aln.data.show, by = list("Group"), FUN = "median")

# ggplot(aln.data.show) +
#   facet_wrap( ~ opstype, scales = "free_x") +
#    geom_boxplot(aes(y = RPM, x = Group, color = median(RPM) > 0), varwidth = T) +
#    #geom_jitter(aes(y = RPM, x = Group), size = .4, width = .3) +
#    #geom_beeswarm(aes(y = RPM, x = Group), size = .1) + 
#    coord_flip() + 
#    theme_light(base_size = 16) #+
# #  stat_summary(fun=sum, colour="darkred", geom="point", shape=18, size=3,show.legend =FALSE)
# #  scale_colour_manual(name = 'median > 0', values = setNames(c('red','green'),c(T, F)))
 
##todo: READS OR READ PAIRS?

# library(dplyr)
# aln.data.g <- group_by(aln.data, Group)
# aln.data.g.s <- summarize(aln.data.g, MWS.RPM = median(MWS.RPM), LWS1.RPM  = median(LWS1.RPM), LWS2.RPM = median(LWS2.RPM))#, nreads = sum(M.Seqs))
# aln.data.gsm <- melt(aln.data.g.s, variable.name = "opstype", value.name = "RPM")

# ggplot(aln.data.show) +
#    facet_wrap( ~ opstype, scales = "free_x") +
#    geom_boxplot(aes(y = RPM, x = Group, color = ifelse(median(aln.data.gsm$RPM) > 0, "black", "red")), varwidth = T) +
#    #geom_jitter(aes(y = RPM, x = Group), size = .4, width = .3) +
#    #geom_beeswarm(aes(y = RPM, x = Group), size = .1) + 
#    coord_flip() + 
#    theme_light(base_size = 16) #+
# 
# 
# 
# mtcars%>%
#    group_by(cyl)%>%
#    mutate(mean_cyl=mean(mpg))#%>%
# #   ggplot(aes(factor(cyl),mpg))+geom_boxplot(aes(fill=mean_cyl))+scale_fill_gradient(low="blue",high = "red")

library(dplyr)

aln.data.show2 <- aln.data.show %>% group_by(Group, opstype) %>% mutate(median = median(RPM)) 

aln.data.show2$Group <- factor(aln.data.show2$Group, levels = c("European", "Baikal 1", "E. cyaneus", "E. verrucosus", "Other Baikal 2",
                                                                "G. lacustris", "G. minus"))


p1 <- 
ggplot(aln.data.show2) +
   facet_wrap( ~ opstype, scales = "fixed") +
   geom_boxplot(aes(y = RPM, x = Group, color = median > 0), varwidth = T) +
   scale_color_manual(values = c("black", "darkviolet"), guide = F) +
   #geom_jitter(aes(y = RPM, x = Group), size = .4, width = .3) +
   #geom_beeswarm(aes(y = RPM, x = Group), size = .1) + 
   coord_flip() + 
   theme_light(base_size = 16) +
   scale_x_discrete(limits = rev(levels(aln.data.show2$Group))) +
   scale_y_continuous(breaks = breaks_pretty(n = 3), limits = c(0, 200)) + 
   theme(axis.text.y=element_text(face = "italic"))
p1


##aln.data.sum <- aln.data %>% group_by(Group) %>% mutate(sum = sum(M.Seqs)) 
aln.data.sum <- aln.data %>% group_by(Group) %>% summarise(sum = sum(M.Seqs)) 
aln.data.sum$dummyfacet <- "Total reads"
##sanity check 
#sum(aln.data[aln.data$Group == "E. cyaneus", "M.Seqs"])


p2 <- ggplot(aln.data.sum) +
   geom_bar(aes(x=Group, y=sum/1000*2), stat="identity") +
   coord_flip() + 
   ylab('Billion reads') + xlab('') +
   theme_light(base_size = 16) +
   facet_wrap(~ dummyfacet) +
   scale_x_discrete(limits = rev(levels(aln.data.show2$Group))) +
   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
##better re-check it
p2   


sum(aln.data[aln.data$Group == "Other Baikal 2", "M.Seqs"])/1000*2
sum(aln.data[aln.data$Group == "Baikal 1", "M.Seqs"])/1000*2

table(aln.data[aln.data$Group == "G. lacustris", "MWS"] == 0)

library(ggpubr)
ggarrange(p1, p2, widths = c(1, .3))

ggsave("alignment_results_scaled.svg", width = 8, height = 4)
#ggsave("alignment_results.png", width = 8, height = 4)


##

 p3 <- 
    ggplot(aln.data.show2) +
    facet_wrap( ~ opstype, scales = "free_x") +
    geom_boxplot(aes(y = RPM, x = Group, color = median > 0), varwidth = T) +
    scale_color_manual(values = c("black", "darkviolet"), guide = F) +
    #geom_jitter(aes(y = RPM, x = Group), size = .4, width = .3) +
    #geom_beeswarm(aes(y = RPM, x = Group), size = .1) + 
    coord_flip() + 
    theme_light(base_size = 14) +
    scale_x_discrete(limits = rev(levels(aln.data.show2$Group))) +
    scale_y_continuous(breaks = breaks_pretty(n = 8)) + 
    theme(axis.text.y=element_text(face = "italic"))
 p3
# 
# ggsave("alignment_results_inserts.svg", width = 6, height = 3)


## inserts one-by-one!

insMWS <- 
   ggplot(aln.data.show2[aln.data.show2$opstype == "MWS.RPM",]) +
#   facet_wrap( ~ opstype, scales = "free_x") +
   geom_boxplot(aes(y = RPM, x = Group, color = median > 0), varwidth = T) +
   scale_color_manual(values = c("black", "darkviolet"), guide = F) +
   #geom_jitter(aes(y = RPM, x = Group), size = .4, width = .3) +
   #geom_beeswarm(aes(y = RPM, x = Group), size = .1) + 
   coord_flip() + 
   theme_light(base_size = 14) +
   scale_x_discrete(limits = rev(levels(aln.data.show2$Group))) +
   scale_y_continuous(breaks = breaks_pretty(n = 1)) + 
   theme(axis.text.y=element_blank(), axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         axis.ticks.y = element_blank())
insMWS

insLWS1 <- 
      ggplot(aln.data.show2[aln.data.show2$opstype == "LWS1.RPM",]) +
      geom_point(aes(y = RPM, x = Group, color = ifelse(median > 0, "darkviolet", "black")))  + 
#      facet_wrap( ~ opstype, scales = "free_x") +
#      geom_boxplot(aes(y = RPM, x = Group, color = ifelse(median > 0, "darkviolet", "black")), varwidth = T) +
      scale_color_manual(values = c("darkviolet", "black"), guide = F) +
      #geom_jitter(aes(y = RPM, x = Group), size = .4, width = .3) +
      #geom_beeswarm(aes(y = RPM, x = Group), size = .1) + 
      coord_flip() + 
      theme_light(base_size = 14) +
      scale_x_discrete(limits = rev(levels(aln.data.show2$Group))) +
      scale_y_continuous(breaks = breaks_pretty(n = 2), limits = c(200, 1200)) + 
      theme(axis.text.y=element_blank(), axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.ticks.y = element_blank())
insLWS1   


insLWS2 <- 
   ggplot(aln.data.show2[aln.data.show2$opstype == "LWS2.RPM",]) +
   geom_point(aes(y = RPM, x = Group, color = ifelse(median > 0, "darkviolet", "black")))  + 
   #      facet_wrap( ~ opstype, scales = "free_x") +
   #      geom_boxplot(aes(y = RPM, x = Group, color = ifelse(median > 0, "darkviolet", "black")), varwidth = T) +
   scale_color_manual(values = c("darkviolet", "black"), guide = F) +
   #geom_jitter(aes(y = RPM, x = Group), size = .4, width = .3) +
   #geom_beeswarm(aes(y = RPM, x = Group), size = .1) + 
   coord_flip() + 
   theme_light(base_size = 14) +
   scale_x_discrete(limits = rev(levels(aln.data.show2$Group))) +
   scale_y_continuous(breaks = breaks_pretty(n = 2), limits = c(200, 600)) + 
   theme(axis.text.y=element_blank(), axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         axis.ticks.y = element_blank())
insLWS2   



ggarrange(p1, p2, widths = c(1, .3))
ggsave("alignment_results_scaled.svg", width = 8, height = 4)
ggarrange(insMWS, insLWS1, insLWS2, nrow = 1) #, width = c(.3, .3, .3))
ggsave("alignment_results_inserts.svg", width = 3, height = 3.3)

