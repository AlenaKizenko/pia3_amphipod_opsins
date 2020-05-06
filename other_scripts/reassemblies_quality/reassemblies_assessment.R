library(tidyr)
library(ggplot2)
df = read.csv('data/busco.csv', stringsAsFactors = F)
df1 = df %>% gather(key = Assembly,value = Busco,-Species)

ggplot(df1, aes(Assembly, Busco)) +
  geom_violin(adjust = 0.8, aes(fill = Assembly)) +
  geom_jitter(alpha=0.5, size=1) +
  theme_classic() +
  geom_boxplot(adjust = 0.7,width=0.15)+
  labs(fill = 'Trinity version', 
        title = 'Comparison of complete BUSCOs percentage in two assemblies',
        y = 'Complete BUSCOs, %',
        x = ' ')+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.title=element_text(size=15), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=15)
        ) +
  scale_x_discrete(labels=c("Complete_BUSCOs_new" = " ", "Complete_BUSCOs_old" = " ")) +
  scale_fill_manual(name = "Trinity version", labels = c('v 2.8.5 (new)', 'v 20140717 (old)'), values = c('seagreen3', 'salmon2'))

old_number = read.csv('data/summary_stats_ass.csv', stringsAsFactors = F)
new_number = read.csv('data/summary_stats_reass.csv', stringsAsFactors = F)

ggplot() +
  geom_histogram(data = old_number, aes(Length, fill = 'Old assemblies'), binwidth = 20, alpha = 0.5) +
  geom_histogram(data = new_number, aes(Length, fill = 'New assemblies'), binwidth = 20, alpha = 0.5) +
  theme_classic()


df_count = read.csv('data/opsins_count.csv', stringsAsFactors = F)
df1_count = df_count %>% gather(key = Assembly,value = Opsins,-Species)
ggplot(df1_count, aes(Assembly, Opsins)) +
  geom_boxplot(width=0.5, aes(fill = Assembly)) +
  theme_classic() +
  labs(fill = 'Trinity version', 
       title = 'Comparison of opsins\' amount in two assemblies',
       y = 'Number of identified opsins',
       x = ' ') +
  theme(legend.position = "none",
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
       # legend.title=element_text(size=15), 
      #  legend.text=element_text(size=12),
        #axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 16),
        axis.title=element_text(size=20)
  ) +
  scale_x_discrete(labels=c("Opsins_number_new" = " ", "Opsins_number_old" = " ")) +
  scale_fill_manual(values = c('seagreen3', 'salmon2'))

op_length = read.csv('data/opsins_length.csv', stringsAsFactors = F)
df1_length = op_length %>% gather(key = Assembly,value = Opsins,-Species)
ggplot(df1_length, aes(Assembly, Opsins)) +
  geom_boxplot(width=0.5, aes(fill = Assembly)) +
  theme_classic() +
  labs(fill = 'Trinity version', 
       title = 'Comparison of opsins\' mean length in two assemblies',
       y = 'Mean length of identified opsins',
       x = ' ') +
  theme(legend.position = "none",
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    #    legend.title=element_text(size=15), 
       # legend.text=element_text(size=12),
       # axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 16),
        axis.title=element_text(size=20)
  ) +
  scale_x_discrete(labels=c("Mean_length_new" = " ", "Mean_length_old" = " ")) +
  scale_fill_manual(values = c('seagreen3', 'salmon2'))
