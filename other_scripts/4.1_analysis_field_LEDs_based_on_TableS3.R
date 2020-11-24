library(openxlsx)

## read the data (Table S3)
results <- read.xlsx("Table_S3_field_exp_LEDs.xlsx")


allInDiods <- sum(results[results$Diod.color %in% c("B", "G", "Y", "R", "K"), "Quantity"])

inBlack <- sum(results[results$Diod.color == "K", "Quantity"])
inBlue <- sum(results[results$Diod.color == "B", "Quantity"])
inGreen <- sum(results[results$Diod.color == "G", "Quantity"])
inYellow <- sum(results[results$Diod.color == "Y", "Quantity"])
inRed <- sum(results[results$Diod.color == "R", "Quantity"])

## sum(results[results$Diod.color == "hand net", "Quantity"])
## 276 in hand net

inBlack / allInDiods * 100 ## 0.37% in black
inBlue / allInDiods * 100 ## 19.1%
inGreen / allInDiods * 100 ## 25.2%
inYellow / allInDiods * 100 ## 22.9%
inRed / allInDiods * 100 ## 32.5%

