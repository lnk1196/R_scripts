library(tidyverse)
library(fs)
library(readxl)
library(ggplot2)
library(data.table)


setwd("~/Desktop")

metadata <- read_excel("PF_metadata-2.xlsx")
completeness <- read_excel("Completeness - Phylofisher.xlsx")

comb <- merge(metadata, completeness, by = "Unique ID")
comb <- subset(comb, select = -c(11:14))
comb <- subset(comb, select = -c(13:15))
setnames(comb, old = c("Unique ID", "Long Name.x", "Higher Taxonomy.x", "Lower Taxonomy.x", "Data Type", "Source", "Raw Data", "Assemblies", "CDS", "Proteomes", "Completeness", "Include in Subset"), 
         new = c("id", "name", "higher_tax", "lower_tax", "data", "source", "raw", "assemblies", "CDS", "proteome", "completeness", "include"))

all_trans <- comb[comb$data == "Transcriptomic",]
trans_CDS <- all_trans[all_trans$CDS == 'x', ]
trans_CDS <- trans_CDS %>% drop_na(CDS)

write.table(trans_CDS, "transciptome_w_CDS.csv", sep = ",")
