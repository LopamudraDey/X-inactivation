

library(dplyr)
library(ggplot2)
library(Matrix)






xci <- read.delim("5227117xci_summary.tsv", sep = "\t")
# Basic cleanup
xci <- xci[nchar(xci$cell) > 3, ]
xci <- xci[xci$cell != "-" & !is.na(xci$cell), ]
xci <- xci[!duplicated(xci$cell), ]


high_conf_xci <- xci[xci$total_snps >= 15, ]

#  Plot both to compare side-by-side
par(mfrow = c(1, 2))

# Original Plot
hist(xci$XCI_score, breaks = 50, 
     main = "Raw XCI (All Cells)", 
     xlab = "XCI score", col = "lightgray")

# Filtered Plot
hist(high_conf_xci$XCI_score, breaks = 50, 
     main = "Filtered XCI (Total SNPs >= 15)", 
     xlab = "XCI score", col = "skyblue")

# 1. Check data structure first
print("Data summary:")
print(summary(xci$total_snps))

# 2. Check how many cells at different thresholds
print(paste("Total cells:", nrow(xci)))
print(paste("Cells with >= 5 SNPs:", sum(xci$total_snps >= 5)))
print(paste("Cells with >= 8 SNPs:", sum(xci$total_snps >= 8)))
print(paste("Cells with >= 10 SNPs:", sum(xci$total_snps >= 10)))

# 3. Choose a realistic threshold  
threshold <- 8 

high_conf_xci <- xci[xci$total_snps >= threshold, ]

# 4. Safely plot only if we actually have cells left!
if (nrow(high_conf_xci) > 0) {
  hist(high_conf_xci$XCI_score, breaks = 30, 
       main = paste("Filtered XCI (Total SNPs >=", threshold, ")"), 
       xlab = "XCI score", col = "skyblue")
} else {
  print("Error: Still no cells left. Your coverage per cell is too low for this threshold.")
}
