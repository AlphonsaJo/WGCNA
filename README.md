# WGCNA

### Weighted Gene Co-Expression Network Analysis (WGCNA)

**Definition**:  
Weighted Gene Co-Expression Network Analysis (WGCNA) is a systems biology method for analyzing gene expression data. It identifies clusters (modules) of highly correlated genes and relates these modules to external traits, helping uncover biological insights, identify biomarkers, and understand functional pathways.

---

## Project Steps for WGCNA

```R
# Step 1: Install and Load Required Libraries
install.packages("WGCNA")
install.packages("ggplot2")
install.packages("dplyr")

library(WGCNA)
library(ggplot2)
library(dplyr)

# Step 2: Prepare and Load Dataset
# Load the dataset
gene_data <- read.csv("dummy_gene_expression_data.csv", row.names = 1)

# Transpose the dataset for WGCNA
datExpr <- t(gene_data)

# Handle missing data
datExpr <- na.omit(datExpr)

# Step 3: Choose Soft-Thresholding Power
powers <- c(1:20)  # Test a range of powers
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot scale independence to choose power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     type = "n", xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit",
     main = "Scale Independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")

# Choose the power value where the model fit reaches ~0.8
softPower <- 6  # Example chosen power, replace based on the plot

# Step 4: Construct Adjacency and TOM Matrix
adjacency <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Step 5: Perform Hierarchical Clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Gene Clustering on TOM-based Dissimilarity", sub = "", xlab = "")

# Step 6: Identify Modules
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE)
dynamicColors <- labels2colors(dynamicMods)

# Plot dendrogram with module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Step 7: Relate Modules to External Traits (Optional)
# Example: Create a dummy trait dataset
traitData <- data.frame(Sample = rownames(datExpr),
                        Trait = c(rep("Control", 5), rep("Disease", 5)))

# Match samples and calculate correlations
moduleTraitCor <- cor(t(datExpr), as.numeric(factor(traitData$Trait)))

# Step 8: Visualize Module-Trait Relationships
heatmap(cor(moduleTraitCor), Rowv = NA, Colv = NA, scale = "column",
        main = "Module-Trait Relationships")
