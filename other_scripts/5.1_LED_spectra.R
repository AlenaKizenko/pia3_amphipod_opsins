library(ggplot2)
library(scales)

B <- read.delim("data/5.1_blue.txt", skip=16, head=F, dec=",")
G <- read.delim("data/5.1_green.txt", skip=16, head=F, dec=",")
Y <- read.delim("data/5.1_yellow.txt", skip=16, head=F, dec=",")
R <- read.delim("data/5.1_red.txt", skip=16, head=F, dec=",")

## Plot the wavelengths
ggplot(data=B, aes(x=V1, y=V2)) + 
  ylab ("Intensity, a.u.") + xlab ("Wavelength, nm") +
  geom_line(data = B, col = "#00C2F9") +
  geom_line(data = G, col = "#00D302") +
  geom_line(data = Y, col = "#FFAC38") +
  geom_line(data = R, col = "#CD022D") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_vline(xintercept = c(380, 760), linetype = "dotted") +  ## approx. borders of the human-visible light spectrum
  theme_bw(base_size = 14)
ggsave("4 spectra.png", width = 9, height = 3)
ggsave("4 spectra.svg", width = 9, height = 3)

## Bandwidths 
## FWHM = full width at half maximal
## https://www.rp-photonics.com/bandwidth.html
fwhm <- function(spectrum) {
  halfintensity <- max(spectrum$V2)/2
  minvalue <- min(spectrum[spectrum$V2 > halfintensity, "V1"])
  maxvalue <- max(spectrum[spectrum$V2 > halfintensity, "V1"])
  return(c(halfintensity, minvalue, maxvalue, maxvalue-minvalue))
}

## Ranges for all LEDs
fwhm(B); fwhm(G); fwhm(Y); fwhm(R)

## Add center wavelength?
centerB <- B$V1[which.max(B$V2)]
centerG <- G$V1[which.max(G$V2)]
centerY <- G$V1[which.max(Y$V2)]
centerR <- G$V1[which.max(R$V2)]

ggplot(data=data.frame(V1=NA, V2=NA), aes(x=V1, y=V2)) + 
  ylab ("Intensity, counts") + xlab ("Wavelength, nm") +
  ## spectrum for the blue LED
  geom_line(data = B, col = "#00C2F9") +
  ## full width at half maximum for the blue LED
  geom_segment(x = fwhm(B)[2], xend = fwhm(B)[3], y=fwhm(B)[1], yend=fwhm(B)[1], col = "#00C2F9") + 
  geom_segment(x = centerB, xend = centerB, y = -10000, yend = max(B$V2), col = "#00C2F9", linetype = "dotted") +
  geom_line(data = G, col = "#00D302") +
  geom_segment(inherit.aes = FALSE, aes(x = fwhm(G)[2], xend = fwhm(G)[3], y=fwhm(G)[1], yend=fwhm(G)[1]), col = "#00D302") + 
  geom_segment(x = centerG, xend = centerG, y = -10000, yend = max(G$V2), col = "#00D302", linetype = "dotted") +
  geom_line(data = Y, col = "#FFAC38") +
  geom_segment(inherit.aes = FALSE, aes(x = fwhm(Y)[2], xend = fwhm(Y)[3], y=fwhm(Y)[1], yend=fwhm(Y)[1]), col = "#FFAC38") + 
  geom_segment(x = centerY, xend = centerY, y = -10000, yend = max(Y$V2), col = "#FFAC38", linetype = "dotted") +
  geom_line(data = R, col = "#CD022D") +
  geom_segment(inherit.aes = FALSE, aes(x = fwhm(R)[2], xend = fwhm(R)[3], y=fwhm(R)[1], yend=fwhm(R)[1]), col = "#CD022D") + 
  geom_segment(x = centerR, xend = centerR, y = -10000, yend = max(R$V2), col = "#CD022D", linetype = "dotted") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
#  scale_x_continuous(breaks = c(380, round(centerB), round(centerG), round(centerY), round(centerR), 760)) +
  annotate(geom = "text", x = c(round(centerB), round(centerG), round(centerY), round(centerR)), 
           y = c(32500, 9500, 18000, 26000), 
           label = c(round(centerB), round(centerG), round(centerY), round(centerR))) + 
  geom_vline(xintercept = c(380, 760), linetype = "dotted") +  ## approx. borders of the human-visible light spectrum
  theme_bw(base_size = 14)

ggsave("4 spectra.png", width = 8, height = 3)
ggsave("4 spectra.svg", width = 8, height = 3)
