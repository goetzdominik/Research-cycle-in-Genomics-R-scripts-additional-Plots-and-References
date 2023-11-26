library(ggplot2)
library(stats)
library(dplyr)
pheno <- read.csv("./pheno.csv", header = TRUE, sep = ";")
#load data with important columns
S68_data <- subset(pheno, pheno[,1] == "Strain6850", select = c("evolution_condition", "growth_in_TSB", "growth_in_supernatant", "average_relative_growth_pqs", "average_STX", "average_H2O2_survival", "Hemolysis_score", "strain_ID"))

#clean columns so that they fit with the genotyp data strain_ID
S68_data[, 8] <- gsub("\\.", "_", S68_data[, 8]) # to get the same ID in both data sets
xT <- read.table("./xT.csv") #gets the rownames of genotyp data
row_names_xT <- rownames(xT)
S68_strings <- S68_data[, 8]

# Initialize a vector to store the positions of non-overlapping row names
row_names_xT <- substr(row_names_xT, 7, nchar(row_names_xT)) # remove strain specific ID, we dont need that
mismatch <- character(0)

# Find where we have data for both datasets, since some gneotype data is missing
for (s1 in S68_strings) {
  match <- FALSE
  for (s2 in row_names_xT) {
    if (s1 == s2) {
      # Set the match flag to TRUE
      match <- TRUE
      # Exit the inner loop since a match is found
      break
    }
  }
  
  # If no match is found, append to mismatch vector
  if (!match) {
    mismatch <- c(mismatch, s1)
  }
}

rows_to_keep <- !(S68_data$strain_ID %in% mismatch)

# Create a new data frame S68_clean containing only the rows to keep
S68_clean <- S68_data[rows_to_keep, ] # created data frame to only have important info for genotyp data
order_row_names <- order(match(S68_clean$strain_ID, row_names_xT))

# Reorder the data frame based on the order of row_names_xT
#created new data frame to make sure that we match the right entries
S68_clean_ordered <- S68_clean[order_row_names, , drop = FALSE]
S68_data <- S68_clean_ordered

# PCA only on important variables 
S68_conti <- subset(S68_data, select = c("growth_in_supernatant" , "growth_in_TSB", "average_relative_growth_pqs", "average_STX", "average_H2O2_survival"))

# Ensure that the selected columns are numeric
S68_conti[, 1:5] <- lapply(S68_conti[, 1:5], function(x) as.numeric(gsub("[^0-9.]", ".", x)))
S68_conti[3:5] <- log2(S68_conti[3:5] + 1)  #log data, because only realtive values given

# Check for and handle missing values
S68_conti <- na.omit(S68_conti)

# Perform PCA on continuous variables
pca_S68_conti <- prcomp(S68_conti, scale = TRUE, center = TRUE)# test if center changes
summary(pca_S68_conti)

#### PCA


# Create a PCA plot for "Strain6850" with color based on "evolution_condition"
#ggplot(data = as.data.frame(pca_S68_conti$x), aes(x = PC1, y = PC2, color = S68_data$evolution_condition)) +
# geom_point() +
#ggtitle("PCA Plot for Strain6850")

#merged data in order to color 
merged_data <- cbind(pca_S68_conti$x,S68_data)
set.seed(123)
# cluster phenotyp data,into 3 clusters
tobe_clustered <- merged_data[,1:2]
cluster2 <- hclust(dist(tobe_clustered))
pheno_cluster <- cutree(cluster2, k = 3)

merged_data <- mutate(merged_data, pheno_cluster)


## originial, without jitter
ggplot(data = merged_data, aes(x = PC1, y = PC2,color = as.character(pheno_cluster), shape = evolution_condition)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "green")) +
  ggtitle("PCA Plot for Strain6850, pheno cluster2")
### PCA PLOT
ggplot(data = merged_data, aes(x = PC1, y = PC2,
                               color = as.character(pheno_cluster), shape = evolution_condition)) +
  geom_point(size = 3, stroke = 0.5, alpha = 0.9, position = position_jitter(width = 0.05, height = 0.05)) +
  scale_color_manual(values = c("#4393c3", "#e74c3c", "#2ecc71"),
                     name = "Phenotype Clusters",
                     labels = c("High grower", "Low grower", "High grower")) +
  scale_shape_manual(name = "Evolution Condition",
                     labels = c("PA Supernatant", "TSB"),
                     values = c(8, 17)) +
  labs(title = "PCA Plot using Phenotype Data",
       x = "PC1",
       y = "PC2") +
  theme_linedraw(base_size = 14) +
  theme(legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 14))


# Assuming your_data is your data frame export cluster
selected_columns <- merged_data[, c("strain_ID", "pheno_cluster")]
write.csv(selected_columns, file = "cluster.csv", row.names = FALSE)

#### LOADINGS FOR THE PCs

#PC2
loadings_S68 <- pca_S68_conti$rotation[, 2]
loadings_data_S68 <- data.frame(Variable = colnames(S68_conti), Loading = loadings_S68)
loadings_data_S68 <- loadings_data_S68[order(abs(loadings_data_S68$Loading), decreasing = TRUE), ]

#plot
ggplot(data = loadings_data_S68, aes(x = Variable, y = Loading, fill = Variable)) +
  geom_bar(stat = "identity", color = "black", fill = "#265783") +
  coord_flip() +
  labs(title = "Variable Loadings for PC2 (Strain6850)",
       x = "yasdsdasd") +  # Removed y-axis label
  theme_linedraw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_blank(),  # Removed y-axis label
        plot.title = element_text(size = 14))



# PC1
loadings_S68_1 <- pca_S68_conti$rotation[,1]

loadings_data_S68_1 <- data.frame(Variable = colnames(S68_conti), Loading = loadings_S68_1)
loadings_data_S68_1 <- loadings_data_S68_1[order(abs(loadings_data_S68_1$Loading), decreasing = TRUE), ]


# plto
ggplot(data = loadings_data_S68_1, aes(x = Variable, y = Loading, fill = Variable)) +
  geom_bar(stat = "identity", color = "black", fill = "#265783") +
  coord_flip() +
  labs(title = "Variable Loadings for PC1 (Strain6850)",
       x = "Loading Type") +  # Removed y-axis label
  theme_linedraw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_blank(),  # Removed y-axis label
        plot.title = element_text(size = 14))








