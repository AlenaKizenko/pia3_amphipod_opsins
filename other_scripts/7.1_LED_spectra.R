library(ggplot2)

B <- read.delim("data/7.1_blue.txt", skip=16, head=F, dec=",")
G <- read.delim("data/7.1_green.txt", skip=16, head=F, dec=",")
Y <- read.delim("data/7.1_yellow.txt", skip=16, head=F, dec=",")
R <- read.delim("data/7.1_red.txt", skip=16, head=F, dec=",")

## Plot the wavelengths
ggplot(data=B, aes(x=V1, y=V2)) + 
  ylab ("Intensity, a.u.") + xlab ("Wavelength, nm") +
  geom_line(data = B, col = "#00C2F9") +
  geom_line(data = G, col = "#00D302") +
  geom_line(data = Y, col = "#FFAC38") +
  geom_line(data = R, col = "#CD022D") +
  geom_vline(xintercept = c(380, 760), linetype = "dotted") +  ## approx. borders of the human-visible light spectrum
  theme_bw(base_size = 14)
ggsave("4 spectra.png", width = 8, height = 3)
ggsave("4 spectra.svg", width = 8, height = 3)

