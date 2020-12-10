library(ggplot2)

## sam files obtained as in 3.1 (bowtie2 & samtools)
## coverage extracted with UGENE

## G. minus (Fig. S3C)

covDat <- read.delim("data/3.3_Gminus_MWScov.txt", header = F)
ggplot(covDat[covDat$V1 == "MWS_Gammarus_pulex_m_NODE_19291_length_2134_cov_5.310647_g14406_i0", ], 
       aes(x=V2, y=V3)) + theme_classic() +
  geom_area(fill="darkgrey") + geom_line(col = "darkgrey") + 
  xlab("Coordinate") + ylab("Coverage") + ggtitle("Alignment of G. minus reads to the G. pulex MWS opsin")
ggsave("Gmi_to_Gpu_MWS.svg", width = 6, height = 4)


## LWS coverage for comparison
# ggplot(covDat[covDat$V1 == "LWS_Gammarus_pulex_m_NODE_27828_length_1449_cov_397.851312_g21126_i0", ], 
#        aes(x=V2, y=V3)) + theme_classic() +
#   geom_area(fill="grey") + geom_line(col = "grey") +
#   xlab("Coordinate") + ylab("Coverage") + ggtitle("Alignment of G. minus reads to a G. pulex LWS opsin")


## G. lacustris (Fig. S3D)

library(ggplot2)

covDat <- read.delim("data/3.3_Glacustris_MWScov.txt", header = F)
ggplot(covDat[covDat$V1 == "MWS_Gammarus_pulex_m_NODE_19291_length_2134_cov_5.310647_g14406_i0", ], 
       aes(x=V2, y=V3)) + theme_classic() +
  geom_area(fill="darkgrey") + geom_line(col = "darkgrey") + 
  xlab("Coordinate") + ylab("Coverage") + ggtitle("Alignment of G. lacustris reads to the G. pulex MWS opsin")
ggsave("Gla_to_Gpu_MWS.svg", width = 6, height = 4)


## LWS coverage for comparison
# ggplot(covDat[covDat$V1 == "LWS_Gammarus_pulex_m_NODE_27828_length_1449_cov_397.851312_g21126_i0", ], 
#        aes(x=V2, y=V3)) + theme_classic() +
#   geom_area(fill="grey") + geom_line(col = "grey") +
#   xlab("Coordinate") + ylab("Coverage") + ggtitle("Alignment of G. minus reads to a G. pulex LWS opsin")



## H. gigas (Fig. S4B)

covDat <- read.delim("data/4.3_Hgigas_hind_LWScov.txt", header = F)
ggplot(covDat, aes(x=V2, y=V3)) + theme_classic(base_size = 14) +
  geom_area(fill="darkgrey") + geom_line(col = "darkgrey") + 
  xlab("Coordinate") + ylab("Coverage") + ggtitle("Alignment of H. gigas reads to the H. gigas LWS opsin")
ggsave("Hgigas_hind_to_Hgigas.svg", width = 8, height = 3)


