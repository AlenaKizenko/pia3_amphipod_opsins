library(ggplot2)
library(ggbeeswarm)
library(openxlsx)

## data input
ovsd <- read.xlsx("busco_both_with_depth.xlsx")
ovsd_nonzero <- ovsd[ovsd$N_LWS_opsins_merged > 0, ]

summary(glm(data = ovsd_nonzero, formula = N_LWS_opsins_merged ~ Typical_depth + Complete_BUSCOs))
#summary(glm(data = ovsd, formula = N_LWS_opsins_merged ~ Typical_depth + Complete_BUSCOs))
summary(glm(data = ovsd_nonzero, formula = N_LWS_opsins_merged ~ log(Typical_depth) + Complete_BUSCOs))
#summary(glm(data = ovsd_nonzero, formula = N_LWS_opsins_merged ~ Typical_mindepth + Complete_BUSCOs + Group))

summary(lm(data = ovsd_nonzero, formula = N_LWS_opsins_merged ~ Typical_depth + Complete_BUSCOs))
summary(lm(data = ovsd_nonzero, formula = N_LWS_opsins_merged ~ log(Typical_depth) + Complete_BUSCOs))


summary(lm(data = ovsd_nonzero, formula = N_LWS_opsins_merged ~ Typical_depth))

hist(ovsd$Complete_BUSCOs, breaks = 60)
hist(ovsd_nonzero$Complete_BUSCOs, breaks = 60)

ggplot(data = ovsd_nonzero, aes(x = N_LWS_opsins_merged, y = Typical_depth, shape = Group)) + 
  geom_beeswarm(alpha = .5, size = 1.5, groupOnX = T, cex = 2, fill = "#6929C4") + 
  theme_bw(base_size = 14) + 
  #scale_shape_manual(values = c(3, 5)) +
  scale_shape_manual(values = c(24, 21)) +
  ylim(1000, 0) + 
  ylab('Typical average depth') + xlab('Number of found LWS opsins') + 
  theme(legend.position = c(0.85, 0.2))

ggsave("Depth_vs_opsins.png", width=12, height=8, units = "cm")
ggsave("Depth_vs_opsins.svg", width=12, height=8, units = "cm")


## BUSCOs and why nopsins = 0 was excluded
ggplot(data = ovsd, aes(x = N_LWS_opsins_merged, y = Complete_BUSCOs, shape = Group)) + 
  geom_beeswarm(alpha = .5, size = 1.5, groupOnX = T, cex = 2) + 
  theme_bw(base_size = 12) + 
  scale_shape_manual(values = c(3, 5)) +
  ylab('Complete BUSCOs') + xlab('Number of found LWS opsins') + 
  theme(legend.position = c(0.85, 0.2))
