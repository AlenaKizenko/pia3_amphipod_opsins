library(reshape2)
library(scales)
library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(ggpubr)

## data input (it's Table S2)
aln.data <- read.csv("TableS2_alignment_to_Gpu_Gammaridae.csv")
## select the data to show
aln.data <- droplevels(aln.data[aln.data$Include.to.main.fig,])

## melt and further rearrange
aln.data.show <- melt(aln.data[, c("Group", "M.Seqs", "MWS.RPM", "LWS1.RPM", "LWS2.RPM")],
                      variable.name = "opstype", value.name = "RPM", id.vars = c("Group", "M.Seqs"))
aln.data.show2 <- aln.data.show %>% group_by(Group, opstype) %>% mutate(median = median(RPM)) 
aln.data.show2$Group <- factor(aln.data.show2$Group, levels = c("European", "Baikal 1", "E. cyaneus", "E. verrucosus", "Other Baikal 2",
                                                                "G. lacustris", "G. minus"))

## Now visualize the alignment results in one scale
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


## The number of reads
aln.data.sum <- aln.data %>% group_by(Group) %>% summarise(sum = sum(M.Seqs)) 
aln.data.sum$dummyfacet <- "Total reads"

p2 <- ggplot(aln.data.sum) +
   geom_bar(aes(x=Group, y=sum/1000*2), stat="identity") +
   coord_flip() + 
   ylab('Billion reads') + xlab('') +
   theme_light(base_size = 16) +
   facet_wrap(~ dummyfacet) +
   scale_x_discrete(limits = rev(levels(aln.data.show2$Group))) +
   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p2   



## Now the inserts, one-by-one
insMWS <- 
   ggplot(aln.data.show2[aln.data.show2$opstype == "MWS.RPM",]) +
   geom_boxplot(aes(y = RPM, x = Group, color = median > 0), varwidth = T) +
   scale_color_manual(values = c("black", "darkviolet"), guide = F) +
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
      scale_color_manual(values = c("darkviolet", "black"), guide = F) +
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
   scale_color_manual(values = c("darkviolet", "black"), guide = F) +
   coord_flip() + 
   theme_light(base_size = 14) +
   scale_x_discrete(limits = rev(levels(aln.data.show2$Group))) +
   scale_y_continuous(breaks = breaks_pretty(n = 2), limits = c(200, 600)) + 
   theme(axis.text.y=element_blank(), axis.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         axis.ticks.y = element_blank())
insLWS2   


## Save the main part and the inserts
ggarrange(p1, p2, widths = c(1, .3))
ggsave("alignment_results_scaled.svg", width = 8, height = 4)
ggarrange(insMWS, insLWS1, insLWS2, nrow = 1) #, width = c(.3, .3, .3))
ggsave("alignment_results_inserts.svg", width = 3, height = 3.3)

## The figures were finally combined in Inkscape
