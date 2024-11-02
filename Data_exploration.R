if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")
BiocManager::install("oligo")
BiocManager::install("limma")
BiocManager::install("GEOquery")
BiocManager::install("tidyverse")

########################################
#  RMA normalization
########################################

# load the libraries

library(affy)
library(GEOquery)
library(tidyverse)

###############################
# for help with the analysis
# openVignette( "affy" )
# library( limma )
# limmaUsersGuide()
###############################

#set the working directory
setwd("C:/Users/evapa/Desktop/AppBio/ΔΙΠΛΩΜΑΤΙΚΗ/Data/RAW")

# get supplementary files for both experiments
getGEOSuppFiles("GSE85543")
getGEOSuppFiles("GSE85544")

# untar files
untar("00Raw_Data/GSE85543/GSE85543_RAW.tar", exdir = '01Untar_data/GSE85543')
untar("00Raw_Data/GSE85544/GSE85544_RAW.tar", exdir = '01Untar_data/GSE85544')

# reading in .cel files
# dir<-"../01Untar_data"
read.data3 <- ReadAffy(celfile.path = "01Untar_data/GSE85543")
read.data4 <- ReadAffy(celfile.path = "01Untar_data/GSE85544")

# performing RMA normalization
normalized.data3 <- rma(read.data3)
normalized.data4 <- rma(read.data4)

# get expression estimates
normalized.expr3 <- as.data.frame(exprs(normalized.data3))
normalized.expr4 <- as.data.frame(exprs(normalized.data4))

# map probe IDs to gene symbols
gse3 <- getGEO("GSE85543", GSEMatrix = TRUE)
gse4 <- getGEO("GSE85544", GSEMatrix = TRUE)

# fetch feature data to get ID - gene symbol mapping
feature.data3 <- gse3$GSE85543_series_matrix.txt.gz@featureData@data
feature.data4 <- gse4$GSE85544_series_matrix.txt.gz@featureData@data
# subset
feature.data3 <- feature.data3[,c(1,11)]
feature.data4 <- feature.data4[,c(1,11)]

#MERGE THE EXPRESSION DATAFRAMES WITH THE FEATURE GENE SYMBOLS TO ALSO HAVE THE EXPRESSION VALUES

normalized.expr3 <- normalized.expr3 %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature.data3, by = 'ID')
# FIRST, SEE IF YOU CAN MOVE IN HERE, THE GENE COLUMN IN THE FRONT
# Move the last column to the place of the first column
normalized.expr3 <- normalized.expr3[, c(ncol(normalized.expr3), 1:(ncol(normalized.expr3) - 1))]
# SECOND REMOVE PROBE COLUMN
normalized.expr3 <- normalized.expr3[-c(2)]

normalized.expr4 <- normalized.expr4 %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature.data4, by = 'ID')
# FIRST, SEE IF YOU CAN MOVE IN HERE, THE GENE COLUMN IN THE FRONT
# Move the last column to the place of the first column
normalized.expr4 <- normalized.expr4[, c(ncol(normalized.expr4), 1:(ncol(normalized.expr4) - 1))]
# SECOND REMOVE PROBE COLUMN
normalized.expr4 <- normalized.expr4[-c(2)]


#TURN IT INTO A TXT FILE
write.table(normalized.expr3, file = "02Normalised_Data/GSE85543_Normalised_expression_data.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(normalized.expr4, file = "02Normalised_Data/GSE85544_Normalised_expression_data.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)


###############################
# PLOT TIME !!!
###############################

###############################
# MA Plots
###############################

# Assuming `exprs_data` is your normalized expression matrix (genes in rows, samples in columns)

# library(limma)

# Load necessary library for plotting
library(ggplot2)

# Exclude the first column (gene names) to get only numeric columns for calculations
numeric_data3 <- normalized.expr3[, -1]

# Calculate the average expression across samples (A value) and the log fold-change (M value)
A_values <- rowMeans(numeric_data3, na.rm = TRUE)
M_values <- apply(numeric_data3, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

# Combine A and M values into a data frame
ma_data3 <- data.frame(A = A_values, M = M_values)

# Create the MA plot
ggplot(ma_data3, aes(x = A, y = M)) +
  geom_point(alpha = 0.5, color = "blue") +
  labs(title = "MA Plot", x = "Average Intensity (A)", y = "Log Fold Change (M)") +
  theme_minimal()

# data4

# Exclude the first column (gene names) to get only numeric columns for calculations
numeric_data4 <- normalized.expr4[, -1]

# Calculate the average expression across samples (A value) and the log fold-change (M value)
A_values <- rowMeans(numeric_data4, na.rm = TRUE)
M_values <- apply(numeric_data4, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

# Combine A and M values into a data frame
ma_data4 <- data.frame(A = A_values, M = M_values)

# Create the MA plot
ggplot(ma_data4, aes(x = A, y = M)) +
  geom_point(alpha = 0.5, color = "blue") +
  labs(title = "MA Plot", x = "Average Intensity (A)", y = "Log Fold Change (M)") +
  theme_minimal()


###############################
# Boxplot
###############################

boxplot(numeric_data3, main = "Boxplot of Normalized Expression Data 3",
        xlab = "Samples", ylab = "Expression Levels",
        col = "lightblue", outline = FALSE)

boxplot(numeric_data4, main = "Boxplot of Normalized Expression Data 4",
        xlab = "Samples", ylab = "Expression Levels",
        col = "lightblue", outline = FALSE) 

#############################
# HISTOGRAM
#############################

# Histogram of M values
ggplot(ma_data3, aes(x = M)) +
  geom_histogram(binwidth = 0.5, fill = "blue", alpha = 0.5, color = "black") +
  labs(title = "Distribution of Log Fold Changes (M)", x = "Log Fold Change (M)", y = "Frequency") +
  theme_minimal()

# Histogram of M values
ggplot(ma_data4, aes(x = M)) +
  geom_histogram(binwidth = 0.5, fill = "blue", alpha = 0.5, color = "black") +
  labs(title = "Distribution of Log Fold Changes (M)", x = "Log Fold Change (M)", y = "Frequency") +
  theme_minimal()

###############################
# Density Plot
###############################

# Load necessary libraries
library(ggplot2)
library(reshape2)

# Converting data to long format for ggplot2
exprs_long3 <- reshape2::melt(normalized.expr3)

ggplot(exprs_long3, aes(value, color = variable)) +
  geom_density() +
  labs(title = "Density Plot of Normalized Expression Data 3 WITH GENES",
       x = "Expression Level", y = "Density") +
  theme_minimal()

# Convert to long format
exprs_long4 <- reshape2::melt(normalized.expr4)

# Create the density plot
ggplot(exprs_long4, aes(value, color = variable)) +
  geom_density() +
  labs(title = "Density Plot of Normalized Expression Data 4 WITH GENES",
       x = "Expression Level", y = "Density") +
  theme_minimal()


##########################
# WITHOUT THE GENE COLUMN

# Converting data to long format for ggplot2
exprs_long3 <- reshape2::melt(numeric_data3)

ggplot(exprs_long3, aes(value, color = variable)) +
  geom_density() +
  labs(title = "Density Plot of Normalized Expression Data 3",
       x = "Expression Level", y = "Density") +
  theme_minimal()

# Convert to long format
exprs_long4 <- reshape2::melt(numeric_data4)

# Create the density plot
ggplot(exprs_long4, aes(value, color = variable)) +
  geom_density() +
  labs(title = "Density Plot of Normalized Expression Data 4",
       x = "Expression Level", y = "Density") +
  theme_minimal()


###############################
# QUALITY CONTROL
###############################

###############################
# Remove NA values
###############################

# FIND THE TOTAL NA VALUES IN THE WHOLE DATAFRAME

total_na3 <- sum(is.na(normalized.expr3))
cat("Total number of NA values:", total_na3, "\n")

total_na4 <- sum(is.na(normalized.expr4))
cat("Total number of NA values:", total_na4, "\n")

# says total number of na values equals 0


# Remove rows with any NA values ?
# clean_data <- na.omit(exprs_data)

# Replace NAs with the median value for each row ?
# exprs_data[is.na(exprs_data)] <- apply(exprs_data, 1, function(x) median(x, na.rm = TRUE))

# Check if there are still any NAs ?
# sum(is.na(exprs_data))

###############################
# Filter out low-expression genes
###############################

# Set the minimum expression threshold (e.g., log2(10))
min_threshold <- log2(10)

# Filter out low-expression genes
filtered.expr3a <- normalized.expr3a[rowMeans(normalized.expr3a[, -ncol(normalized.expr3a)], na.rm = TRUE) > min_threshold, ]
filtered.expr4a <- normalized.expr4a[rowMeans(normalized.expr4a[, -ncol(normalized.expr4a)], na.rm = TRUE) > min_threshold, ]

# Check the number of genes before and after filtering
cat("Number of genes before filtering (GSE85543):", nrow(normalized.expr3a), "\n")
cat("Number of genes after filtering (GSE85543):", nrow(filtered.expr3a), "\n")
cat("Number of genes before filtering (GSE85544):", nrow(normalized.expr4a), "\n")
cat("Number of genes after filtering (GSE85544):", nrow(filtered.expr4a), "\n")

###############################
# Outlier Detection
###############################

# Calculate Z-scores
z_scores3a <- scale(filtered.expr3a[, -ncol(filtered.expr3a)])  # Exclude probeset ID column
z_scores4a <- scale(filtered.expr4a[, -ncol(filtered.expr4a)])  # Exclude probeset ID column

# Set Z-score threshold for outlier detection (e.g., > 3 or < -3)
outlier_threshold <- 3

# Identify outliers
outliers3a <- which(abs(z_scores3a) > outlier_threshold, arr.ind = TRUE)
outliers4a <- which(abs(z_scores4a) > outlier_threshold, arr.ind = TRUE)

# Extract indices of outlier genes
outlier_genes3a <- unique(rownames(filtered.expr3a)[outliers3a[, 1]])
outlier_genes4a <- unique(rownames(filtered.expr4a)[outliers4a[, 1]])

# Output the number of outliers detected
cat("Number of outlier genes detected (GSE85543):", length(outlier_genes3a), "\n")
cat("Number of outlier genes detected (GSE85544):", length(outlier_genes4a), "\n")

###############################
# PLOT TIME !!! AGAIN !!!
###############################

###############################
# MA Plots
###############################

# Assuming `exprs_data` is your normalized expression matrix (genes in rows, samples in columns)

library(limma)

# Example for the first two columns
ma_data <- exprs_data[, c(1, 2)] # Select two columns for comparison
plotMA(ma_data, main = "MA Plot", ylim = c(-5, 5))

###############################
# Boxplot
###############################

boxplot(exprs_data, main = "Boxplot of Normalized Expression Data",
        xlab = "Samples", ylab = "Expression Levels",
        col = "lightblue", outline = FALSE)

###############################
# Density Plot
###############################

library(ggplot2)

# Converting data to long format for ggplot2
exprs_long <- reshape2::melt(exprs_data)

ggplot(exprs_long, aes(value, color = Var2)) +
  geom_density() +
  labs(title = "Density Plot of Normalized Expression Data",
       x = "Expression Level", y = "Density") +
  theme_minimal()

###############################
# RLE and NUSE plots, if available -> check what these are
###############################

###############################
# NEXT STEPS: PCA & HEATMAP & VOLCANO
###############################
