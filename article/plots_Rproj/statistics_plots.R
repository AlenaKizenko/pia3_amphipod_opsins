library(ggplot2)
library(ggpubr)
library(wesanderson)
library(psycho)
library(cowplot)

#------------------------reading file with statistics--------------------------------------------------------------------------
busco = read.csv('data/busco_both.csv', stringsAsFactors = F)

#------------testing data normality-----------------------------------------------------------------------------------------
shapiro.test(busco$Single_BUSCOs)
shapiro.test(busco$Single_BUSCOs)
shapiro.test(busco$Single_BUSCOs)
shapiro.test(busco$Single_BUSCOs)
shapiro.test(busco$Single_BUSCOs)
shapiro.test(busco$Single_BUSCOs)


#--------------------building Complete BUSCOs percentage plot----------------------------------------------------------------------------
complete_b = ggplot(busco, aes(Assembler, Complete_BUSCOs, fill = Assembler)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(alpha=0.7, size=2, position = position_jitter(height = .05, width = .1)) +
  theme_bw() +
  labs(fill = 'Assembler', 
       y = 'Percentage',
       x = '') +
  theme(plot.title = element_text(size = 20),
        legend.position = "none",
        legend.title=element_text(size=16), 
        legend.text=element_text(size=16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20)) +
  scale_fill_manual(values = wes_palette("IsleofDogs1", n = 2)) +
 stat_compare_means(method = "t.test", label.y = 100, size = 4) +
 ggtitle("Complete BUSCOs")

#--------------------building Single copy BUSCOs percentage plot----------------------------------------------------------------------------
single_b = ggplot(busco, aes(Assembler, Single_BUSCOs, fill = Assembler)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(alpha=0.7, size=2, position = position_jitter(height = .05, width = .1)) +
  theme_bw() +
  labs(fill = 'Assembler', 
       y = '',
       x = '') +
  theme(plot.title = element_text(size = 20),
        legend.position = "right",
        legend.title=element_text(size=16), 
        legend.text=element_text(size=16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20)) +
  scale_fill_manual(values = wes_palette("IsleofDogs1", n = 2)) +
  stat_compare_means(method = "t.test", label.y = 90, size = 4) +
  ggtitle("Single copy BUSCOs")

#--------------------building Duplicated BUSCOs percentage plot----------------------------------------------------------------------------
duplicated_b = ggplot(busco, aes(Assembler, Duplicated_BUSCOs, fill = Assembler)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(alpha=0.7, size=2, position = position_jitter(height = .05, width = .1)) +
  theme_bw() +
  labs(fill = 'Assembler', 
       y = 'Percentage',
       x = '') +
  theme(plot.title = element_text(size = 20),
        legend.position = "none",
        legend.title=element_text(size=16), 
        legend.text=element_text(size=16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20)) +
  scale_fill_manual(values = wes_palette("IsleofDogs1", n = 2)) +
  stat_compare_means(method = "t.test", label.y = 60, size = 4) +
  ggtitle("Duplicated BUSCOs")

#--------------------building Fragmented BUSCOs percentage plot----------------------------------------------------------------------------
fragmented_b = duplicated_b = ggplot(busco, aes(Assembler, Fragmented_BUSCOs, fill = Assembler)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(alpha=0.7, size=2, position = position_jitter(height = .05, width = .1)) +
  theme_bw() +
  labs(fill = 'Assembler', 
       y = 'Percentage',
       x = '') +
  theme(plot.title = element_text(size = 20),
        legend.position = "none",
        legend.title=element_text(size=16), 
        legend.text=element_text(size=16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20)) +
  scale_fill_manual(values = wes_palette("IsleofDogs1", n = 2)) +
  stat_compare_means(method = "t.test", label.y = 35, size = 4) +
  ggtitle("Fragmented BUSCOs")

#--------------------building Missing BUSCOs percentage plot--------------------------------------------------------------------------
missing_b = duplicated_b = ggplot(busco, aes(Assembler, Missing_BUSCOs, fill = Assembler)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(alpha=0.7, size=2, position = position_jitter(height = .05, width = .1)) +
  theme_bw() +
  labs(fill = 'Assembler', 
       y = '',
       x = '') +
  theme(plot.title = element_text(size = 20),
        legend.position = "right",
        legend.title=element_text(size=16), 
        legend.text=element_text(size=16),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20)) +
  scale_fill_manual(values = wes_palette("IsleofDogs1", n = 2)) +
  stat_compare_means(method = "t.test", label.y = 65, size = 4) +
  ggtitle("Missing BUSCOs")

#------------------------combining Complete, Single copy and Duplicated BUSCOs in one plot--------------------------------------------
plot_grid(complete_b + scale_y_continuous(limits = c(0, 100), n.breaks = 5),
          single_b + scale_y_continuous(limits = c(0, 100), n.breaks = 5),
          duplicated_b + scale_y_continuous(limits = c(0, 100), n.breaks = 5),
          rel_heights = c(1,1),
          rel_widths = c(1,1.3))

#------------------------combining Fragmented and Missing BUSCOs in one plot--------------------------------------------
plot_grid(fragmented_b + scale_y_continuous(limits = c(0, 100), n.breaks = 5),
          missing_b + scale_y_continuous(limits = c(0, 100), n.breaks = 5),
          rel_heights = c(1,1),
          rel_widths = c(1,1.4))

