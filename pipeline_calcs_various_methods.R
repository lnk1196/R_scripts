library(tidyverse)
library(fs)
library(readxl)
library(ggplot2)

setwd("~/Desktop/Pipeline Data Files/")

# Create an empty data frame to store the results
SNPdata <- data.frame(abbrev = character(), SNPden = numeric(), stringsAsFactors = FALSE)
SNP_R <- data.frame(abbrev = character(), SNPden = numeric(), stringsAsFactors = FALSE)
SNPamoeba <- data.frame(abbrev = character(), SNPden = numeric(), stringsAsFactors = FALSE)
SNPcrypto <- data.frame(abbrev = character(), SNPden = numeric(), stringsAsFactors = FALSE)

# Read the Excel file
sheets <- excel_sheets("SNPpipeline_knowns.xlsx")
sheets_restricted <- excel_sheets("Refined_SNP_data.xlsx")
sheets_amoebae <- excel_sheets("SNPpipeline_amoebae.xlsx")
sheets_crypto <- excel_sheets("SNPpipeline_Cyrpto.xlsx")

# Loop through each sheet
for (sheet in sheets) {
  # Read the sheet
  df <- read_excel("SNPpipeline_knowns.xlsx", sheet = sheet)
  
  # Calculate the average of the SNP_COUNT column
  avg_snp <- mean(df$SNP_COUNT, na.rm = TRUE)
  
  # Calculate the number of unique names in the first column
  num_unique_names <- length(unique(df$CHROM))
  
  # Create a new row with sheet name, average, and number of unique names
  new_row <- data.frame(abbrev = sheet, SNPden = avg_snp, Contigs_w_SNP = num_unique_names, stringsAsFactors = FALSE)
  
  # Add the new row to SNPdata data frame
  SNPdata <- rbind(SNPdata, new_row)
}

for (sheet in sheets_restricted) {
  df <- read_excel("Refined_SNP_data.xlsx", sheet = sheet)
  avg_snp <- mean(df$SNP_COUNT, na.rm = TRUE)
  num_unique_names <- length(unique(df$CHROM))
  new_row <- data.frame(abbrev = sheet, SNPden = avg_snp, Contigs_w_SNP = num_unique_names, stringsAsFactors = FALSE)
  SNP_R <- rbind(SNP_R, new_row)
}

for (sheet in sheets_amoebae) {
  df <- read_excel("SNPpipeline_amoebae.xlsx", sheet = sheet)
  avg_snp <- mean(df$SNP_COUNT, na.rm = TRUE)
  num_unique_names <- length(unique(df$CHROM))
  new_row <- data.frame(abbrev = sheet, SNPden = avg_snp, Contigs_w_SNP = num_unique_names, stringsAsFactors = FALSE)
  SNPamoeba <- rbind(SNPamoeba, new_row)
}

for (sheet in sheets_crypto) {
  df <- read_excel("SNPpipeline_Cyrpto.xlsx", sheet = sheet)
  avg_snp <- mean(df$SNP_COUNT, na.rm = TRUE)
  num_unique_names <- length(unique(df$CHROM))
  new_row <- data.frame(abbrev = sheet, SNPden = avg_snp, Contigs_w_SNP = num_unique_names, stringsAsFactors = FALSE)
  SNPcrypto <- rbind(SNPcrypto, new_row)
}

# Remove unneeded rows
SNPdata <- SNPdata[-(1:3),]
SNP_R <- SNP_R[-(1:3),]
SNPamoeba <- SNPamoeba[-(1:3),]
SNPcrypto <- SNPcrypto[-(1:3),]

# Import Analysis sheet 
dat_df <- read_excel("SNPpipeline_knowns.xlsx", 2)
rest_df <- read_excel("Refined_SNP_data.xlsx", 2)
amoe_df <- read_excel("SNPpipeline_amoebae.xlsx", 2)
cryp_df <- read_excel("SNPpipeline_Cyrpto.xlsx", 2)

# Copy columns to SNPdata
SNPdata$Ploidy <- dat_df$Ploidy
SNPdata$SNPs <- dat_df$SNPs
SNPdata$Bases <- dat_df$Bases
SNPdata$Total_Contigs <- dat_df$Total_Contigs
SNPdata$Type <- rep("F", nrow(SNPdata))

# Copy columns to SNP_R
SNP_R$Ploidy <- rest_df$Ploidy
SNP_R$SNPs <- rest_df$SNPs
SNP_R$Bases <- rest_df$Bases
SNP_R$Total_Contigs <- rest_df$Total_Contigs
SNP_R$Type <- rep("R", nrow(SNP_R))

# Copy columns to SNPamoeba
SNPamoeba$Ploidy <- amoe_df$Ploidy
SNPamoeba$SNPs <- amoe_df$SNPs
SNPamoeba$Bases <- amoe_df$Bases
SNPamoeba$Total_Contigs <- amoe_df$Total_Contigs
SNPamoeba$Type <- amoe_df$Type

# Copy columns to SNPcrypto
SNPcrypto$Ploidy <- cryp_df$Ploidy
SNPcrypto$SNPs <- cryp_df$SNPs
SNPcrypto$Bases <- cryp_df$Bases
SNPcrypto$Total_Contigs <- cryp_df$Total_Contigs
SNPcrypto$Type <- cryp_df$Type

# Calculations SNPdata
SNPdata$SNP_KB <- (SNPdata$SNPs / SNPdata$Bases)*1000
SNPdata$wSNP <- (SNPdata$Contigs_w_SNP/SNPdata$Total_Contigs)
SNPdata$woSNP <- (1 - SNPdata$wSNP)

# Calculations SNP_R
SNP_R$SNP_KB <- (SNP_R$SNPs / SNP_R$Bases)*1000
SNP_R$wSNP <- (SNP_R$Contigs_w_SNP/SNP_R$Total_Contigs)
SNP_R$woSNP <- (1 - SNP_R$wSNP)

# Calculations SNPamoeba
SNPamoeba$SNP_KB <- (SNPamoeba$SNPs / SNPamoeba$Bases)*1000
SNPamoeba$wSNP <- (SNPamoeba$Contigs_w_SNP/SNPamoeba$Total_Contigs)
SNPamoeba$woSNP <- (1 - SNPamoeba$wSNP)

# Calculations SNPcrypto
SNPcrypto$SNP_KB <- (SNPcrypto$SNPs / SNPcrypto$Bases)*1000
SNPcrypto$wSNP <- (SNPcrypto$Contigs_w_SNP/SNPcrypto$Total_Contigs)
SNPcrypto$woSNP <- (1 - SNPcrypto$wSNP)

# Save calc values to CVS
write.csv(SNPdata, file = "SNPdenCalc_knowns.csv")
write.csv(SNP_R, file = "SNPdenCalc_refined.csv")
write.csv(SNPamoeba, file = "SNPdenCalc_amoebae.csv")
write.csv(SNPcrypto, file = "SNPdencalc_Crypoper.csv")

#Combine into a single df
all_data_df <- bind_rows(SNPdata, SNP_R, SNPamoeba, SNPcrypto)

restricted <- subset(all_data_df, Type == "R")
restricted <- restricted[!duplicated(restricted$abbrev), ]

full_method <- subset(all_data_df, Type == "F")
full_method <- full_method[!duplicated(full_method$abbrev), ]
full_method$PlotOrder <- ifelse(full_method$Ploidy == "unknown", 2, 1)
full_method <- full_method[order(full_method$PlotOrder), ]

known_F <- subset(full_method, Ploidy != "unknown")
known_F <- known_F[!duplicated(known_F$abbrev), ]

write_csv(all_data_df, file = "SNPdenCalc_all.csv")

##########################################################
# Plotting

#Only plot orgs with known ploidy
ggplot(known_F)+geom_bar(aes(x=factor(abbrev, level=c("Scere1", "Scere2", "Amoeprot", "Caraaura", "Solatube", "Tritaest", "Fraganan")), y = woSNP), position = "fill", stat = "identity", width = 0.5) + geom_bar(aes(x = abbrev, y = wSNP, fill = Ploidy), position = "stack", stat = "identity", width = 0.5)  + labs(x = "", y = "Percentage of Contigs with SNPs")
ggplot(known_F)+geom_bar(aes(x=factor(abbrev, level=c("Scere1", "Scere2", "Amoeprot", "Caraaura", "Solatube", "Tritaest", "Fraganan")), y = SNPden, fill = Ploidy), stat = "identity", width = 0.5) + labs(x = "", y = "SNPs per kilobase")
ggplot(data = known_F) + geom_bar(aes(x = factor(abbrev, level=c("Scere1", "Scere2", "Amoeprot", "Caraaura", "Solatube", "Tritaest", "Fraganan")), y = SNP_KB, fill = Ploidy), stat = "identity", width = 0.5) + labs(x = "", y = "Calculated SNPs per kilobase")
ggplot(data = known_F) + geom_bar(aes(x = factor(abbrev, level=c("Scere1", "Scere2", "Amoeprot", "Caraaura", "Solatube", "Tritaest", "Fraganan")), y = Total_Contigs), stat = "identity", width = 0.5) + geom_bar(aes(x = abbrev, y = Contigs_w_SNP, fill = Ploidy), stat = "identity", width = 0.5) + labs(x = "", y = "Ratio of Contigs with SNPs to Total Contigs")

#Only plot orgs with known ploidy from Phylofisher method
ggplot(SNP_R)+geom_bar(aes(x=factor(abbrev, level=c("Scere1", "Scere2", "Amoeprot", "Caraaura", "Solatube", "Tritaest", "Fraganan")), y = SNPden, fill = Ploidy), stat = "identity", width = 0.5) + labs(title = "Pipeline SNP density",x = "", y = "SNPs per kilobase")
ggplot(data = SNP_R) + geom_bar(aes(x = factor(abbrev, level=c("Scere1", "Scere2", "Amoeprot", "Caraaura", "Solatube", "Tritaest", "Fraganan")), y = SNP_KB, fill = Ploidy), stat = "identity", width = 0.5) + labs(title = "Calculated SNP density", x = "", y = "Calculated SNPs per kilobase")
ggplot(data = SNP_R) + geom_bar(aes(x = factor(abbrev, level=c("Scere1", "Scere2", "Amoeprot", "Caraaura", "Solatube", "Tritaest", "Fraganan")), y = Total_Contigs), stat = "identity", width = 0.5) + geom_bar(aes(x = abbrev, y = Contigs_w_SNP, fill = Ploidy), stat = "identity", width = 0.5) + labs(x = "", y = "Contigs with SNPs vs. Total Contigs")
ggplot(SNP_R)+geom_bar(aes(x=factor(abbrev, level=c("Scere1", "Scere2", "Amoeprot", "Caraaura", "Solatube", "Tritaest", "Fraganan")), y = woSNP),position = "fill", stat = "identity", width = 0.5) + geom_bar(aes(x = abbrev, y = wSNP, fill = Ploidy), position = "stack", stat = "identity", width = 0.5)  + labs(title = "Percent of Total Contigs containing SNPS", x = "", y = "Percentage of Contigs with SNPs")

#Only plot data for the amoebae
ggplot(SNPamoeba) + geom_bar(aes(x = abbrev, y = SNPden, fill = Ploidy), stat = "identity", width = 0.5) + labs(x = "", y = "SNPs per kilobase")
ggplot(data = SNPamoeba) + geom_bar(aes(x = abbrev, y = SNP_KB, fill = Ploidy), stat = "identity", width = 0.5) + labs(x = "", y = "Calculated SNPs per kilobase")
ggplot(data = SNPamoeba) + geom_bar(aes(x = abbrev, y = Total_Contigs), stat = "identity", width = 0.5) + geom_bar(aes(x = abbrev, y = Contigs_w_SNP, fill = Ploidy), stat = "identity", width = 0.5) + labs(x = "", y = "Ratio of Contigs with SNPs to Total Contigs")
ggplot(SNPamoeba)+geom_bar(aes(x=abbrev, y = woSNP),position = "fill", stat = "identity", width = 0.5) + geom_bar(aes(x = abbrev, y = wSNP, fill = Ploidy), position = "stack", stat = "identity", width = 0.5)  + labs(x = "", y = "Percentage of Contigs with SNPs")

# Plot data for crypto
ggplot(SNPcrypto) + geom_bar(aes(x = abbrev, y = SNPden, fill = Type), stat = "identity", width = 0.5) + labs(x = "", y = "SNPs per kilobase") + facet_wrap(~Type, scales = "free_y")
ggplot(data = SNPcrypto) + geom_bar(aes(x = abbrev, y = SNP_KB, fill = Type), stat = "identity", width = 0.5) + labs(x = "", y = "Calculated SNPs per kilobase") + facet_wrap(~Type, scales = "free_y")
ggplot(data = SNPcrypto) + geom_bar(aes(x = abbrev, y = Total_Contigs), stat = "identity", width = 0.5) + geom_bar(aes(x = abbrev, y = Contigs_w_SNP, fill = Type), stat = "identity", width = 0.5) + labs(x = "", y = "Ratio of Contigs with SNPs to Total Contigs") + facet_wrap(~Type, scales = "free_y")
ggplot(SNPcrypto)+geom_bar(aes(x=abbrev, y = woSNP),position = "fill", stat = "identity", width = 0.5) + geom_bar(aes(x = abbrev, y = wSNP, fill = Type), position = "stack", stat = "identity", width = 0.5)  + labs(x = "", y = "Percentage of Contigs with SNPs") + facet_wrap(~Type, scales = "free_x")

# Full Transcriptome Method
ggplot(full_method)+geom_bar(aes(x= abbrev, y = woSNP), position = "fill", stat = "identity", width = 0.5) + geom_bar(aes(x = abbrev, y = wSNP, fill = Ploidy), position = "stack", stat = "identity", width = 0.5)  + labs(x = "", y = "Percentage of Contigs with SNPs")
ggplot(full_method)+geom_bar(aes(x= abbrev, y = SNPden, fill = Ploidy), stat = "identity", width = 0.5) + labs(x = "", y = "SNPs per kilobase")
ggplot(full_method) + geom_bar(aes(x = abbrev, y = SNP_KB, fill = Ploidy), stat = "identity", width = 0.5) + labs(x = "", y = "Calculated SNPs per kilobase")
ggplot(data = full_method) + geom_bar(aes(x = abbrev, y = Total_Contigs), stat = "identity", width = 0.5) + geom_bar(aes(x = abbrev, y = Contigs_w_SNP, fill = Ploidy), stat = "identity", width = 0.5) + labs(x = "", y = "Ratio of Contigs with SNPs to Total Contigs")


# All organisms
ggplot(full_method, aes(x = reorder(abbrev, PlotOrder), y = SNPden, fill = Ploidy)) + geom_bar(stat = "identity") + labs(title = "SNP Density by Organism", x = "Organism", y = "SNP Density")
ggplot(full_method, aes(x = reorder(abbrev, PlotOrder), y = SNP_KB, fill = Ploidy)) + geom_bar(stat = "identity") + labs(title = "Calculated SNP Density by Organism", x = "Organism", y = "SNP Density")
ggplot(data = full_method) + geom_bar(aes(x = reorder(abbrev, PlotOrder), y = woSNP), position = "fill", stat = "identity", width = 0.5) + geom_bar(aes(x = reorder(abbrev, PlotOrder), y = wSNP, fill = Ploidy), position = "stack", stat = "identity", width = 0.5) + labs(title = "Percent of Total Contigs with SNPS", x = "Organism", y = "Contig Percentage")


#Plot knowns together and unknowns together side by side
known_ploidy_data <- full_method[full_method$Ploidy != "unknown", ]
unknown_ploidy_data <- full_method[full_method$Ploidy == "unknown", ]
ggplot() + geom_bar(data = known_ploidy_data, aes(x=factor(abbrev, level=c("Scere1", "Scere2", "Amoeprot", "Caraaura", "Solatube", "Tritaest", "Fraganan")), y = SNPden, fill = Ploidy), stat = "identity", width = 0.5) + labs(title = "SNP Density by Organism", y = "SNP Density") + facet_wrap(~PlotOrder, scales = "free_x") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_bar(data = unknown_ploidy_data, aes(x = reorder(abbrev, PlotOrder), y = SNPden, fill = Ploidy), stat = "identity", position = "identity", width = 0.5) + labs(x = "Organism")
ggplot() + geom_bar(data = known_ploidy_data, aes(x=factor(abbrev, level=c("Scere1", "Scere2", "Amoeprot", "Caraaura", "Solatube", "Tritaest", "Fraganan")), y = SNP_KB, fill = Ploidy), stat = "identity", width = 0.5) + labs(title = "SNP Density by Organism", y = "SNP Density") + facet_wrap(~PlotOrder, scales = "free_x") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_bar(data = unknown_ploidy_data, aes(x = reorder(abbrev, PlotOrder), y = SNP_KB, fill = Ploidy), stat = "identity", position = "identity", width = 0.5) + labs(x = "Organism")

ggplot() + geom_bar(data = known_ploidy_data, aes(x=factor(abbrev, level=c("Scere1", "Scere2", "Amoeprot", "Caraaura", "Solatube", "Tritaest", "Fraganan")), y = Contigs_w_SNP, fill = Ploidy), stat = "identity", width = 0.5) + geom_bar(aes(x = abbrev, y = Total_Contigs), position = "stack", stat = "identity", width = 0.5) + labs(title = "SNP Density by Organism", y = "SNP Density") + facet_wrap(~PlotOrder, scales = "free_x") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_bar(data = unknown_ploidy_data, aes(x = reorder(abbrev, PlotOrder), y = Contigs_w_SNP, fill = Ploidy), stat = "identity", position = "identity", width = 0.5) + geom_bar(aes(x = reorder(abbrev, PlotOrder), y = Total_Contigs), position = "stack", stat = "identity", width = 0.5)+labs(x = "Organism")


## Faceted plots for Total Contigs vs Contigs w SNPs

# Define the specific order for known_ploidy_data
specific_order <- c("Scere1", "Scere2", "Amoeprot", "Caraaura", "Solatube", "Tritaest", "Fraganan")
known_ploidy_data$abbrev <- factor(known_ploidy_data$abbrev, levels = specific_order)

# Combine the data frames
combined_data <- rbind(known_ploidy_data, unknown_ploidy_data)

# Plotting using ggplot2 with facet_wrap
ggplot(combined_data) +
  geom_bar(aes(x = reorder(abbrev, PlotOrder), y = Total_Contigs), 
           position = "stack", stat = "identity", width = 0.5) +
  geom_bar(aes(x = reorder(abbrev, PlotOrder), y = Contigs_w_SNP, fill = Ploidy), 
           position = "stack", stat = "identity", width = 0.5) +
  labs(title = "Number of Contigs Containing SNPs", y = "Contigs", x = "Organism") +
  facet_wrap(~PlotOrder, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(combined_data) +
  geom_bar(aes(x = reorder(abbrev, PlotOrder), y = woSNP), 
           position = "fill", stat = "identity", width = 0.5) +
  geom_bar(aes(x = reorder(abbrev, PlotOrder), y = wSNP, fill = Ploidy), 
           position = "stack", stat = "identity", width = 0.5) +
  labs(title = "Percentage of Contigs Containing SNPs", y = "Percent", x = "Organism") +
  facet_wrap(~PlotOrder, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(combined_data) +
  geom_bar(aes(x = reorder(abbrev, PlotOrder), y = Total_Contigs), 
           position = "fill", stat = "identity", width = 0.5) +
  geom_bar(aes(x = reorder(abbrev, PlotOrder), y = percentage, fill = Ploidy), 
           position = "stack", stat = "identity", width = 0.5) +
  labs(title = "Percentage of Contigs Containing SNPs", y = "Percent", x = "Organism") +
  facet_wrap(~PlotOrder, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

