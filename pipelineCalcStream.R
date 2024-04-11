# Load libraries
library(tidyverse)
library(readxl)

# Set working directory
setwd("~/Desktop/Pipeline Data Files/")

# Read data
read_data <- function(file_path, sheet_index) {
  sheet_name <- excel_sheets(file_path)[sheet_index]
  read_excel(file_path, sheet = sheet_name)
}

# Data processing
process_data <- function(data) {
  data$SNPs <- as.numeric(data$SNPs)
  data$Bases <- as.numeric(data$Bases)
  data$Contigs_w_SNP <- as.numeric(data$Contigs_w_SNP)
  data$Total_Contigs <- as.numeric(data$Total_Contigs)

  data$SNP_KB <- (data$SNPs / data$Bases) * 1000
  data$wSNP <- (data$Contigs_w_SNP / data$Total_Contigs)
  data$woSNP <- (1 - data$wSNP)
  
  data
}

# Save data to CSV
save_to_csv <- function(data, filename) {
  write.csv(data, file = filename)
}

# Plotting functions
plot_SNP_density <- function(data) {
  ggplot(data, aes(x = Abbrev, y = SNPden, fill = Ploidy)) +
    geom_bar(stat = "identity", width = 0.5) +
    labs(x = "", y = "SNPs per kilobase") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(labels = scales::number_format(scale = 0.01))
}

# Define other plotting functions similarly...

# Save plots to files
save_plots <- function(plot, filename) {
  ggsave(filename, plot, width = 10, height = 6)
}

# Main script

# Specify file path and sheet index
file_path <- "SNPpipeline_All.xlsx"
sheet_index <- 2

# Read data
SNP_all <- read_data(file_path, sheet_index)

# Process data
processed_data <- process_data(SNP_all)

# Save processed data to CSV
save_to_csv(processed_data, "SNPdencalc_All.csv")

# Generate and save plots
restricted_SNPden_plot <- plot_SNP_density(processed_data)
save_plots(restricted_SNPden_plot, "restricted_SNPden.png")

# Continue defining other plotting functions and saving plots...

# Filter data for selected samples
selected_samples <- c("Arabthal2", "Arabthal3", "Arabthal4", "Scere1", "Scere2", "Physpate1", "Physpate2", "Physpoly1", "Physpoly2")
restricted_selected <- subset(processed_data, Abbrev %in% selected_samples)

# Generate and save multi-plots
multi_SNPden_plot <- plot_SNP_density(restricted_selected)
save_plots(multi_SNPden_plot, "multi_SNPden.png")

# Continue generating and saving multi-plots...

