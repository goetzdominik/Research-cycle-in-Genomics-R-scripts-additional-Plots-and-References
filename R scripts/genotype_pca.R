library(vcfR)
library(adegenet)
library(dplyr)
library(tidyverse)
library(ggplot2)

vcf <- read.vcfR("./genom data/SA6850/2020-09-25_mod_SA6850_filteredFromAncestor.vcf", verbose = FALSE)
x <- vcfR2genind(vcf, return.alleles = TRUE, ploidy = 1, sep = "/")

xT <- tab(x, freq=TRUE, NA.method="zero")
metadata = as.data.frame(rownames(xT)) %>% 
  mutate(ID = rownames(xT)) %>% 
  separate('rownames(xT)', into = c("clone", "Treatment", "Tube", "Clone")) %>% 
  column_to_rownames("ID")


x = x[rownames(metadata[!metadata$Clone == "P",]),]
#x = x[rownames(metadata[!metadata$Treatment == "T",]),]
xT <- tab(x, freq=TRUE, NA.method="zero")
metadata = as.data.frame(rownames(xT)) %>% 
  mutate(ID = rownames(xT)) %>% 
  separate('rownames(xT)', into = c("clone", "Treatment", "Tube", "Clone")) %>% 
  column_to_rownames("ID")
write.table(xT, file = "xT.csv")
pca1 <- dudi.pca(df = xT, center = TRUE, scale = FALSE, scannf = FALSE, nf = 15)

 #s.label(pca1$li)




cluster <-read.csv("../Blockkurs RCG/cluster.csv")
plot(pca1$eig)
metadata<- cbind(metadata,cluster)

#original
merge(pca1$li, metadata, by = "row.names") %>% 
  ggplot(aes(Axis6, Axis2, color =as.character(pheno_cluster), shape = Treatment)) + geom_point() +
  ggtitle("PCA Plot for Strain6850, Treatment")


merge(pca1$li, metadata, by = "row.names") %>%
  ggplot(aes(Axis6, Axis2, color = as.character(pheno_cluster), shape = Treatment)) +
  geom_point(size = 3, stroke = 0.5, alpha = 0.9, position = position_jitter(width = 0.05, height = 0.05)) +
  scale_color_manual(values =  c("#4393c3", "#e74c3c", "#2ecc71"),
                     name = "Phenotype Clusters",
                     labels = c("High grower", "Low grower", "High grower"))  +
  scale_shape_manual(name = "Evolution Condition",
                     labels = c("PA Supernatant", "TSB"),
                     values = c( 8, 17)) +
  labs(title = "PCA Plot using Genotype Data",
       x = "PC6",
       y = "PC2") +
  theme_linedraw(base_size = 14) +
  theme(legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 13))

library(devtools)
devtools::install_github("riboseq/Ribo")
library(Ribo)

# Set a specific seed for reproducibility
set.seed(123)

# Define the blue color from your previous plots
desired_blue <- "#4393c3"

# Modify loading plot for PC2


loadingplot(pca1$c1, axis = 2, thres = 0.25, lab.jitter = 3.5, main = "PC2",
            col = desired_blue,
            pch = 16,
            cex.lab = 1.5,
            cex.axis = 1.2,
            cex.main = 1.8,
            cex.sub = 1.2,
            xlab="Position in Genome",
            labels.custom = customize_labels,
            custom_names = c("11", "Lab1el2", "Label3"))
            

# Modify loading plot for PC6
loadingplot(pca1$c1, axis = 1, thres = 0.25, lab.jitter = 2.5, main = "PC6",
            col = desired_blue,
            pch = 16,
            cex.lab = 1.5,
            cex.axis = 1.2,
            cex.main = 1.8,
            cex.sub = 1.2,
            xlab="Position in Genome",
            labels.custom = customize_labels,
            custom_names = c("Label4", "Label5", "Label6"))



x <- abs(pca1$c1[,1]) > 0.25
rownames(pca1$c1[x,1])



