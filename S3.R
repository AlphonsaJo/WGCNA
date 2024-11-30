# Advanced WGCNA Analysis with Enhanced Visualizations

# Install and load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("WGCNA", "ComplexHeatmap"))
install.packages(c("ggplot2", "tidyverse", "cluster", "pheatmap", 
                   "umap", "plotly", "factoextra", "ggdendro"))

library(WGCNA)
library(ggplot2)
library(tidyverse)
library(cluster)
library(pheatmap)
library(umap)
library(plotly)
library(factoextra)
library(ComplexHeatmap)
library(ggdendro)

# Set working directory
setwd("D:/Documents/Course_files/Sem7/mmd/MMD Project/Final_Dataset")

# Read the dataset
data <- read.csv("TCGA_rpkm.csv", row.names = 1)

# Data preprocessing
# Remove columns with zero variance
data <- data[, apply(data, 2, var) != 0]

# Normalization
data_normalized <- scale(data)

# 1. WGCNA Network Analysis
# Choose soft thresholding powers
powers <- c(1:20)
sft <- pickSoftThreshold(data_normalized, powerVector = powers, verbose = 5)

# Select optimal soft power
soft_power <- sft$fitIndices[which.max(sft$fitIndices[,2]), 1]

# Network construction
adjacency <- adjacency(data_normalized, power = soft_power)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Hierarchical clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Module detection
dynamicMods <- cutreeDynamic(dendro = geneTree, 
                              distM = dissTOM,
                              deepSplit = 2, 
                              pamRespectsDendro = FALSE,
                              minClusterSize = 30)

# 2. Advanced Visualizations

# 2.1 Interactive Soft Threshold Plot
pdf("interactive_soft_threshold.pdf")
plot_ly(x = sft$fitIndices[,1], y = -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        type = 'scatter', mode = 'markers',
        text = paste("Power:", sft$fitIndices[,1]),
        title = "Scale Free Topology Model Fit")
dev.off()

# 2.2 Comprehensive Heatmap with Clustering
pdf("comprehensive_heatmap.pdf")
pheatmap(data_normalized[1:100, 1:50],  # Sample subset for visualization
         scale = "row", 
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Gene Expression Heatmap")
dev.off()

# 2.3 UMAP Visualization
umap_result <- umap(data_normalized)
pdf("umap_visualization.pdf")
plot(umap_result$layout, col = dynamicMods, 
     pch = 16, main = "UMAP Visualization of Gene Expression")
dev.off()

# 2.4 PCA Visualization with Variance Explained
pca_result <- prcomp(data_normalized)
pdf("pca_variance_plot.pdf")
fviz_eig(pca_result, addlabels = TRUE, 
         main = "Variance Explained by Principal Components")
dev.off()

# 2.5 Dendrogram with Module Colors
# Assign module colors
moduleColors <- labels2colors(dynamicMods)
pdf("colored_dendrogram.pdf")
plotDendroAndColors(geneTree, moduleColors, "Module Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# 2.6 K-means Clustering Visualization
kmeans_result <- kmeans(pca_result$x[,1:2], centers = 5)
pdf("enhanced_kmeans.pdf")
fviz_cluster(kmeans_result, data = pca_result$x[,1:2],
             ellipse.type = "convex",
             palette = "Set2",
             ggtheme = theme_minimal())
dev.off()

# 3. Correlation and Network Analysis
# Module eigengenes
MEs <- moduleEigengenes(data_normalized, dynamicMods)$eigengenes

# Inter-module correlation
pdf("module_correlation.pdf")
cor_modules <- cor(MEs)
pheatmap(cor_modules, 
         main = "Inter-module Correlation",
         color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# 4. Network Topology Analysis
# Connectivity analysis
connectivity <- softConnectivity(data_normalized, power = soft_power)
pdf("connectivity_distribution.pdf")
hist(connectivity, breaks = 100, main = "Connectivity Distribution")
dev.off()

# 5. Summarize Results
sink("analysis_summary.txt")
cat("WGCNA Analysis Summary\n")
cat("====================\n")
cat("Soft Power Threshold:", soft_power, "\n")
cat("Number of Modules:", length(unique(dynamicMods)), "\n")
cat("K-means Clustering Centers:\n")
print(kmeans_result$centers)
sink()

# Save workspace
save.image("comprehensive_analysis_results.RData")
