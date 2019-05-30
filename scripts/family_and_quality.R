library(ggplot2)
df = read.csv('/home/alena/Documents/IB/project_opsins/stats.csv')
str(df)
df$Species = as.character(df$Species)
str(df)
df$Species = as.factor(df$Species)
str(df)
ggplot(df, aes(fill = Family, x = Species, y = Missing.BUSCOs)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  facet_grid(. ~ Family, space = 'free_x', scales = 'free_x') +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  scale_fill_manual(values = c("#FF6666", "#999999", "#FF9966", "#00CC99", "#CC3300", "#CC6699", "#336600", "#FFCC66", "#3399CC"))