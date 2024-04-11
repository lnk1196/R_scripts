library(tidyverse)
library(fs)
library(readxl)
library(ggplot2)

setwd("~/Desktop/Pipeline Data Files/")

SNPdata <- data.frame(abbrev = character(), SNPden = numeric(), stringsAsFactors = FALSE)

  sheets <- excel_sheets("known_rerun.xlsx")
for (sheet in sheets) {
  # Read the sheet
  df <- read_excel("known_rerun.xlsx", sheet = sheet)
  
  # Calculate the average of the SNP_COUNT column
  avg_snp <- mean(df$SNP_COUNT, na.rm = TRUE)
  
  # Calculate the number of unique names in the first column
  num_unique_names <- length(unique(df$CHROM))
  
  # Create a new row with sheet name, average, and number of unique names
  new_row <- data.frame(abbrev = sheet, SNPden = avg_snp, Contigs_w_SNP = num_unique_names, stringsAsFactors = FALSE)
  
  # Add the new row to SNPdata data frame
  SNPdata <- rbind(SNPdata, new_row)
}

#SNPdata <- SNPdata[-(1:3),]

# Import Analysis sheet 
dat_df <- read_excel("known_rerun.xlsx", 1)

# Copy columns to SNPdata
SNPdata$SNPs <- dat_df$snps
SNPdata$Bases <- dat_df$bases
SNPdata$Total_Contigs <- dat_df$Total_Contigs

# Make numeric
dat_df$SNPs <- as.numeric(dat_df$SNPs)
dat_df$Bases <- as.numeric(dat_df$Bases)
dat_df$Contigs_w_SNP <- as.numeric(dat_df$Contigs_w_SNP)
dat_df$Total_Contigs <- as.numeric(dat_df$Total_Contigs)


ggplot(dat_df)+geom_bar(aes(x= Abbrev, y = woSNP), position = "fill", stat = "identity", width = 0.5) + geom_bar(aes(x = Abbrev, y = wSNP, fill = Run), position = "stack", stat = "identity", width = 0.5) + labs(x = "Organism", y = 'Percent of Contigs with SNPs', title = "Comparison of Percent of Contigs with SNPs (new vs. old run)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(dat_df)+geom_bar(aes(x= Abbrev, y = SNPden, fill = Run), stat = "identity", width = 0.5) + labs(x = "Organism", y = 'Pipeline SNP density', title = "Comparison of Pipeline SNP Density (new vs. old run)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(dat_df) + geom_bar(aes(x = Abbrev, y = SNP_KB, fill = Run), stat = "identity", width = 0.5) + labs(x = "Organism", y = 'Calculated SNP density (SNPs/KB)', title = "Comparison of Calculated SNP Density (new vs. old run)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(data = dat_df) + geom_bar(aes(x = Abbrev, y = Total_Contigs), stat = "identity", width = 0.5) + geom_bar(aes(x = Abbrev, y = Contigs_w_SNP, fill = Run), stat = "identity", width = 0.5) + labs(x = "Organism", y = 'Number of Contigs with SNPs', title = "Comparison of Number of Contigs with SNPs (new vs. old run)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
