library(ggplot2)
library(ggpubr)
library(wesanderson) #install.packages("wesanderson")
library(psycho)
library(cowplot)
library(svglite)
library(forcats)   ##for fct_reorder

#------------------------reading file with statistics--------------------------------------------------------------------------
busco = read.csv('data/1.8_busco_both.csv', stringsAsFactors = F)

#------------testing data normality-----------------------------------------------------------------------------------------
shapiro.test(busco$Single_BUSCOs)
shapiro.test(busco$Complete_BUSCOs)
shapiro.test(busco$Duplicated_BUSCOs)
shapiro.test(busco$Fragmented_BUSCOs)
shapiro.test(busco$Missing_BUSCOs)

## A convenience list to plot 
plot.details <- list(geom_boxplot(width = 0.5, outlier.alpha = 0),
                     geom_jitter(alpha=0.7, size=1.5, position = position_jitter(height = .05, width = .1)),
                     theme_bw(),
                     labs(fill = 'Assembler', 
                          y = '',
                          x = ''),
                     theme(plot.title = element_text(size = 14), #20
                           legend.position = "right",
                           legend.title=element_text(size=12), #16
                           legend.text=element_text(size=12), #16
                           axis.text.x = element_text(size = 12, angle=70, vjust = 1, hjust=1), #12
                           axis.text.y = element_text(size = 12), #16
                           axis.title=element_text(size=12)), #16
                     scale_fill_manual(values = wes_palette("IsleofDogs1", n = 3),
                                       labels=c("rnaSPAdes", "Trinity (TSA)", "Trinity 2.8.5")),
                     scale_x_discrete(labels=c("rnaSPAdes" = "rnaSPAdes", "Trinity" = "Trinity (TSA)",
                                               "Trinity285" = "Trinity 2.8.5")),
                     scale_y_continuous(limits = c(0, 1300), n.breaks = 5))
##

#--------------------building Complete BUSCOs percentage plot----------------------------------------------------------------------------
my_comparisons <- compare_means(p.adjust.method = "holm", method="wilcox.test", data = busco, 
                                formula = Complete_BUSCOs ~ Assembler) %>% mutate(Assembler = group1)

complete_b = ggplot(busco, aes(x = fct_reorder(Assembler, Complete_BUSCOs*10.66), y = Complete_BUSCOs*10.66, fill = Assembler)) +
  plot.details +
 ggtitle("Complete BUSCOs")+
  stat_pvalue_manual(my_comparisons, label = "p.adj", y.position = 1000, step.increase = 0.1)
complete_b

#--------------------building Single copy BUSCOs percentage plot----------------------------------------------------------------------------
my_comparisons <- compare_means(p.adjust.method = "holm", method="wilcox.test", data = busco, 
                                formula = Single_BUSCOs ~ Assembler) %>% mutate(Assembler = group1)

single_b = ggplot(busco, aes(x = fct_reorder(Assembler, Complete_BUSCOs*10.66), Single_BUSCOs*10.66, fill = Assembler)) +
  plot.details +
  ggtitle("Single copy BUSCOs") +
  stat_pvalue_manual(my_comparisons, label = "p.adj", y.position = 1000, step.increase = 0.1)

#--------------------building Duplicated BUSCOs percentage plot----------------------------------------------------------------------------
my_comparisons <- compare_means(p.adjust.method = "holm", method="wilcox.test", data = busco, 
                                formula = Duplicated_BUSCOs ~ Assembler) %>% mutate(Assembler = group1)

duplicated_b = ggplot(busco, aes(x = fct_reorder(Assembler, Complete_BUSCOs*10.66), Duplicated_BUSCOs*10.66, fill = Assembler)) +
  plot.details +
  ggtitle("Duplicated BUSCOs") +
  stat_pvalue_manual(my_comparisons, label = "p.adj", y.position = 1000, step.increase = 0.1)


#--------------------building Fragmented BUSCOs percentage plot----------------------------------------------------------------------------
my_comparisons <- compare_means(p.adjust.method = "holm", method="wilcox.test", data = busco, 
                                formula = Fragmented_BUSCOs ~ Assembler) %>% mutate(Assembler = group1)

fragmented_b = ggplot(busco, aes(x = fct_reorder(Assembler, Complete_BUSCOs*10.66), Fragmented_BUSCOs*10.66, fill = Assembler)) +
  plot.details +
  stat_pvalue_manual(my_comparisons, label = "p.adj", y.position = 1000, step.increase = 0.1) +
  ggtitle("Fragmented BUSCOs")

#--------------------building Missing BUSCOs percentage plot--------------------------------------------------------------------------

my_comparisons <- compare_means(p.adjust.method = "holm", method="wilcox.test", data = busco, 
                                formula = Missing_BUSCOs ~ Assembler) %>% mutate(Assembler = group1)

my_comparisons$Assembler <- my_comparisons$group1

missing_b <- ggplot(busco, aes(x = fct_reorder(Assembler, Complete_BUSCOs*10.66), Missing_BUSCOs*10.66, fill = Assembler)) +
  plot.details +
  stat_pvalue_manual(my_comparisons, label = "p.adj", y.position = 1000, step.increase = 0.1) +
  ggtitle("Missing BUSCOs")


#-------------------------------Opsin plot----------------------------------------
my_comparisons <- compare_means(p.adjust.method = "holm", method="wilcox.test", data = busco, 
                                formula = Opsins_found ~ Assembler) %>% mutate(Assembler = group1)

opsins_found = 
  ggplot(busco, aes(x = fct_reorder(Assembler, Complete_BUSCOs*10.66), Opsins_found, fill = Assembler)) +
  geom_boxplot(width = 0.5) +
  geom_point() +
  geom_count() +
  #geom_jitter(alpha=0.7, size=2, position = position_jitter(height = .05, width = .1)) +
  theme_bw() +
  geom_line(aes(group = Species), alpha = .5) +
  labs(fill = 'Assembler', 
       y = '',
       x = '') +
  theme(plot.title = element_text(size = 14),
        legend.position = "right",
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(size = 12, angle=70, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=12)) +
  scale_fill_manual(values = wes_palette("IsleofDogs1", n = 3),
                    labels=c("rnaSPAdes", "Trinity (TSA)", "Trinity 2.8.5")) +
  scale_x_discrete(labels=c("rnaSPAdes" = "rnaSPAdes", "Trinity" = "Trinity (TSA)",
                            "Trinity285" = "Trinity 2.8.5")) +
  stat_pvalue_manual(my_comparisons, label = "p.adj", y.position = 5.5, step.increase = 0.1) +
  expand_limits(y=7) +
  ggtitle("Opsins found")


opsins_found
compare_means(data = busco, formula = Opsins_found ~ Assembler)


plot_grid(complete_b +  theme(legend.position = 'none'),
          single_b +  theme(legend.position = 'none'),
          duplicated_b,# + theme(legend.position = 'none'),
          fragmented_b + theme(legend.position = 'none'),
          missing_b + theme(legend.position = 'none'),
          opsins_found,
          nrow = 2,
          labels = "AUTO",
#          rel_heights = c(1,1,1,1,1),
          rel_widths = c(1,1,1.5))

ggsave("Assemblies.pdf", width = 25, height = 25, units = "cm")
ggsave("Assemblies.svg", width = 25, height = 25, units = "cm")
