library(ggplot2)

B <- read.delim("blue.txt", skip=16, head=F, dec=",")
G <- read.delim("green.txt", skip=16, head=F, dec=",")
Y <- read.delim("yellow.txt", skip=16, head=F, dec=",")
R <- read.delim("red.txt", skip=16, head=F, dec=",")
UV <- read.delim("UV.txt", skip = 16, head = F, dec = ",")
#ggplot(UV, aes(x=V1, y=V2)) + geom_line() + theme_bw()


#UVf <- read.delim("UV + IR filter.txt", skip = 16, head = F, dec = ",")
IR <- read.delim("Orient.txt", skip = 16, head = F, dec = ",")
IR935 <- read.delim("Mi.txt", skip = 16, head = F, dec = ",")

ggplot(UV, aes(x=V1, y=V2)) + geom_line(col = "violet") + 
  ylab ("Intensity, a.u.") + xlab ("Wavelength, nm") +
  #geom_line(data = UVf, col = "darkviolet") +
  geom_line(data = B, col = "#00C2F9") +
  geom_line(data = G, col = "#00D302") +
  geom_line(data = Y, col = "#FFAC38") +
  geom_line(data = R, col = "#CD022D") +
  geom_line(data = IR, col = "red4") +
  geom_vline(xintercept = c(380, 760), linetype = "dotted") + 
  theme_bw(base_size = 14)
ggsave("spectra.png")
  

ggplot(UV, aes(x=V1, y=V2)) + geom_line(col = "violet", size=2) + 
  ylab ("Intensity, a.u.") + xlab ("Wavelength, nm") +
  #geom_line(data = UVf, col = "darkviolet") +
  geom_line(data = B, col = "#00C2F9", size = 2) +
  geom_line(data = G, col = "#00D302", size = 2) +
  geom_line(data = Y, col = "#FFAC38", size = 2) +
  geom_line(data = R, col = "#CD022D", size = 2) +
  xlim(330, 700) +
  theme_bw(base_size = 14)
ggsave("spectra_concl.png")
