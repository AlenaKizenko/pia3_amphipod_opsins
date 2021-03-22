library("phangorn")

## main text (Fig. 5A), Baikal gammarids only
## This nexus file was created with Splitstree 4
network <- read.nexus.networx("data/6.1_Fig5A_BaikalGammaridaeLWSntandGmi.best.fas.best.fas.nex")
nn <- 2446 #nvertices #can get eg with grep nvertices

## Vector of color(s)
lvc <- adjustcolor(rep("#6929C4", nn), alpha.f=1)
## Make a vector to store node symbols
lva <- rep("", nn)
lva[network$translate$node] <- c("\U25CF") ## the default value
##correct for the Baikal 1 group
lva[network$translate$node[c(grep("Micruropus", network$translate$label, fixed=TRUE))]] <- c("\U25B2")
lva[network$translate$node[c(grep("Gmelinoides", network$translate$label, fixed=TRUE))]] <- c("\U25B2")
lva[network$translate$node[c(grep("Crypturopus", network$translate$label, fixed=TRUE))]] <- c("\U25B2")
lva[network$translate$node[c(grep("Linevichella", network$translate$label, fixed=TRUE))]] <- c("\U25B2")
lva[network$translate$node[c(grep("Asprogammarus", network$translate$label, fixed=TRUE))]] <- c("\U25B2")
lva[network$translate$node[c(grep("Baikalogammarus", network$translate$label, fixed=TRUE))]] <- c("\U25B2")
lva[network$translate$node[c(grep("Macrohectopus", network$translate$label, fixed=TRUE))]] <- c("\U25B2")
## and G. minus
lva[network$translate$node[c(grep("minus", network$translate$label, fixed=TRUE))]] <- c("\U2739")
## If only one opsin, according to Table S1, then the symbols should be empty:
lva[network$translate$node[c(grep("Asprogammarus_rhodophthalmus", network$translate$label, fixed=TRUE))]] <- c("\U25B3")
lva[network$translate$node[c(grep("Micruropus_parvulus", network$translate$label, fixed=TRUE))]] <- c("\U25B3")
## The same for the Baikal 2 group:
lva[network$translate$node[c(grep("Boeckaxelia_carpenterii", network$translate$label, fixed=TRUE))]] <- c("\U25CB")
lva[network$translate$node[c(grep("Carinurus_bicarinatus", network$translate$label, fixed=TRUE))]] <- c("\U25CB")
lva[network$translate$node[c(grep("Eulimnogammarus_ussolzewii", network$translate$label, fixed=TRUE))]] <- c("\U25CB")
lva[network$translate$node[c(grep("Macropereiopus_parvus", network$translate$label, fixed=TRUE))]] <- c("\U25CB")
lva[network$translate$node[c(grep("Ommatogammarus_albinus", network$translate$label, fixed=TRUE))]] <- c("\U25CB")
lva[network$translate$node[c(grep("Ommatogammarus_flavus", network$translate$label, fixed=TRUE))]] <- c("\U25CB")
lva[network$translate$node[c(grep("Pachyschesis_branchialis", network$translate$label, fixed=TRUE))]] <- c("\U25CB")
lva[network$translate$node[c(grep("Parapallasea_borowskii", network$translate$label, fixed=TRUE))]] <- c("\U25CB")
lva[network$translate$node[c(grep("Parapallasea_wosnessenskii", network$translate$label, fixed=TRUE))]] <- c("\U25CB")

## Save the result (Figure 5A)
svg("network.svg", width = 5, height = 5)
plot(network, type="2D", show.tip.label=FALSE, show.node.label=TRUE, node.label=lva, 
     edge.width=0.5, cex=1.5, tip.color=adjustcolor(lvc, alpha.f=0.5) )

legend("topleft", border = FALSE, bty = "n" ,
       legend = c("\U25B2 Baikal 1", "\U25CF Baikal 2", "\U2739 G. minus"),
       col = "#6929C4")
dev.off()



## Supplement; all opsins (excluding non-opsins, thus 144 instead of 146)
## read the network file created with Splitstree4
network <- read.nexus.networx("data/6.1_FigS6_144opsins.best.fas.nex")
nn <- 4511 #nvertices #can get eg with grep nvertices

## Vector of colors
lvc <- adjustcolor(rep("black", nn), alpha.f=1)
## Make a vector to store node symbols
lva <- rep("", nn)
lva[network$translate$node] <- c("\U25CF") #c("\U2727") #c("\U2B24") 
lvc <- adjustcolor(rep("#6929C4", nn), alpha.f=1)
## Adjust for MWS 
lva[network$translate$node[c(grep("MWS", network$translate$label, fixed=TRUE))]] <- c("\U2727")
lva[network$translate$node[c(grep("SWS", network$translate$label, fixed=TRUE))]] <- c("\U2602")
lva[network$translate$node[c(grep("VERL", network$translate$label, fixed=TRUE))]] <- c("\U2605")
##Baikal2
lvc[network$translate$node[c(grep("Eulimnogammarus", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("allasea", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Echiuropus", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Hyallelopsis", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Palicarinus", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Heterogammarus", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Brandtia", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Dorogostaiskia", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Pentagonurus", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Poekilogammarus", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Carinurus", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Brachyuropus", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Ommatogammarus", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Brachyuropus", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Eucarinogammarus", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Pachys", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Hyalellopsis", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Odontogammarus", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("Macropereiopus", network$translate$label, fixed=TRUE))]] <- "#0072B2"
lvc[network$translate$node[c(grep("canth", network$translate$label, fixed=TRUE))]] <- "#0072B2"
#Baikal 1
lvc[network$translate$node[c(grep("Gmelinoides", network$translate$label, fixed=TRUE))]] <- "#5fb2e6"
lvc[network$translate$node[c(grep("Micruropus", network$translate$label, fixed=TRUE))]] <- "#5fb2e6"
lvc[network$translate$node[c(grep("Asprogammarus", network$translate$label, fixed=TRUE))]] <- "#5fb2e6"
lvc[network$translate$node[c(grep("Crypturopus", network$translate$label, fixed=TRUE))]] <- "#5fb2e6"
lvc[network$translate$node[c(grep("Linevichella", network$translate$label, fixed=TRUE))]] <- "#5fb2e6"
#All (other) gammarids
lvc[network$translate$node[c(grep("Gammarus", network$translate$label, fixed=TRUE))]] <- "darkgreen"
lvc[network$translate$node[c(grep("Echinogammarus", network$translate$label, fixed=TRUE))]] <- "darkgreen"
lvc[network$translate$node[c(grep("Marinogammarus", network$translate$label, fixed=TRUE))]] <- "darkgreen"
#Talitridae
lvc[network$translate$node[c(grep("Parhyale", network$translate$label, fixed=TRUE))]] <- '#CC79A7'
lvc[network$translate$node[c(grep("Hyalella", network$translate$label, fixed=TRUE))]] <- '#CC79A7'
lvc[network$translate$node[c(grep("Talitrus", network$translate$label, fixed=TRUE))]] <- '#CC79A7'


## Save the result (Figure S5)
svg("all_opsins_network.svg", width = 10, height = 5)
plot(network, type="2D", show.tip.label=FALSE, show.node.label=TRUE, node.label=lva, 
     edge.width=0.5, cex=1, cex.edge.label = 0.05, tip.color = lvc)
legend("topright", border = FALSE, bty = "n" ,
       legend = c("\U2B24 LWS", "\U2727 MWS", "\U2602 UV/SWS", "\U2605 VERL"),
       col = "#6929C4")
legend("bottomright", border = FALSE, bty = "n" ,
       legend = c("Baikal 2", "Baikal 1", "Other Gammaridae", "Talitridae", "Other amphipods"),
       fill = c("#0072B2", "#5fb2e6", "darkgreen", '#CC79A7', "black"))
dev.off()
