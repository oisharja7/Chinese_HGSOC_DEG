if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(limma)

#setting working directory
getwd()
setwd("D:\\DEG analysis_oisharja\\Chinese HGSOC")

# loading files
files <- list.files(pattern = "*.txt")
files

raw_data <- read.maimages(files, source = "agilent", green.only = TRUE)
rownames(raw_data)

head(raw_data$E) #QC

#Although agilent has less bg noise, QC for bg noise through boxplot
boxplot(raw_data$E,
        main = "Raw Intensity Values",
        xlab = "Arrays",
        ylab = "Intensity value",
        col = "green",
        outline = FALSE)

#second QC
plotDensities(raw_data$E,
              main = "Density Plot of Intensity Values",
              col = rainbow(ncol(raw_data$E))) 

#bg noise reduction 
bg_corrected_data <- backgroundCorrect(raw_data, method = "normexp", offset = 50)

boxplot(bg_corrected_data$E,
        main = "Raw Intensity Values",
        xlab = "Arrays",
        ylab = "Intensity value",
        col = "green",
        outline = FALSE)

plotDensities(bg_corrected_data$E,
              main = "Density Plot of Intensity Values",
              col = rainbow(ncol(bg_corrected_data$E))) # I dont know what happened :')


#normalization
normalized_data <- normalizeBetweenArrays(bg_corrected_data, method = "quantile")
head(normalized_data)
rownames(normalized_data$E)[1:10]

#normalization
normalized_data <- normalizeBetweenArrays(bg_corrected_data, method = "quantile")
head(normalized_data)
rownames(normalized_data$E)[1:10]

#reading and merging annotation file 
annotation <- read.table("GPL16699-15607.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
head(annotation)

#reading and merging annotation file 
#Alternative approach (ITS A SUCCESS!)

#Convert normalized_data$E to dataframe
normalized_df <- as.data.frame(normalized_data)

#assign a row number column to new df
normalized_df$ProbeID <- rownames(normalized_df)

#merge 
annoatated_data <- merge(normalized_df, annotation, by.x = "ProbeID", by.y = "ID", all.x = TRUE)
head(annoatated_data)

# DEG 

# vector of sample labels 
sample_labels <- c("Normal", "Tumor", "Normal", "Tumor", 
                   "Normal", "Tumor", "Normal", "Tumor", 
                   "Normal", "Tumor", "Normal", "Tumor", 
                   "Normal", "Tumor", "Normal", "Tumor", 
                   "Normal", "Tumor", "Normal", "Tumor", 
                   "Normal", "Tumor", "Normal", "Tumor", 
                   "Normal", "Tumor", "Normal", "Tumor", 
                   "Normal", "Tumor") 
sample_labels

# factorizing group
group <- factor(sample_labels, levels = c("Normal", "Tumor"))

# building design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
colnames(design)
design

# defining contrast 
contrast <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)

# fit the linear model
fit <- lmFit(normalized_data, design)

# apply contrast
contrast_fit <- contrasts.fit(fit, contrast)

#empirical bayes analysis
contrast_fit <- eBayes(contrast_fit)

#extract DEG
deg_results <- topTable(contrast_fit, coef = 1, adjust.method = "fdr", number = Inf)
deg_results

#significant ones
significant_genes <- subset(deg_results, adj.P.Val < 0.05 & abs(logFC) > 1)
significant_genes

#final QC
hist(deg_results$adj.P.Val)

hist(deg_results$logFC)

#export result 
write.csv(significant_genes, "significant_genes.csv", row.names = FALSE)

#Visualization 

#Volcano Plot
install.packages("ggplot2")

library(ggplot2)

#Add a column for significance 
deg_results$Significance <- ifelse(deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1, 
                                   "Significant", 
                                   "Not Significant")
colnames(deg_results)
table(deg_results$Significance)

#Plot
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("red", "gray")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme_minimal()


#install.packages("pheatmap")

library(pheatmap)

# top 50 significant DEGs
#op_genes <- deg_results[order(deg_results$adj.P.Val), ][1:50, ]
#xpression_matrix <- normalized_data$E[rownames(normalized_data$E) %in% top_genes$GeneID, ]

# Plot heatmap
#heatmap(expression_matrix,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = data.frame(Group = group),
         main = "Heatmap of Top 50 DEGs")

# TROUBLESHOOT
# Check the structure
str(normalized_data$E)

# Inspect the first few rows and columns
#ead(normalized_data$E)
#ownames(normalized_data$E)

# Preview top_genes$GeneID
#ead(top_genes$ProbeName)


# Filter the expression matrix based on ProbeName
# Ensure 'top_genes' and 'normalized_data$E' are properly aligned
#verlap_probes <- rownames(normalized_data$E)[rownames(normalized_data$E) %in% top_genes$ProbeName]
#xpression_matrix <- normalized_data$E[overlap_probes, ]
# Set rownames to the matching ProbeName from top_genes
# This ensures the heatmap displays the correct labels
#ownames(expression_matrix) <- overlap_probes
# Optional: Scale the data for better visualization
#caled_matrix <- t(scale(t(expression_matrix)))  # Scale rows (probes)

# Create the heatmap
#heatmap(
  expression_matrix               # Scaled expression matrix
  cluster_rows = TRUE,           # Cluster by rows
  cluster_cols = TRUE,           # Cluster by columns
  show_rownames = TRUE,          # Show probe names as row labels
  show_colnames = TRUE,          # Show sample names as column labels
  color = colorRampPalette(c("blue", "white", "red"))(50),  # Color gradient
  fontsize_row = 6,              # Adjust font size for rows
  fontsize_col = 8,              # Adjust font size for columns
  main = "Heatmap of Selected Probes"  # Title of the heatmap
)
#PCA plot

library(ggplot2)

# Perform PCA
pca <- prcomp(t(normalized_data$E), scale. = TRUE)

# Prepare data for plotting
pca_data <- as.data.frame(pca$x)
pca_data$Group <- group

# Plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme_minimal()

head(significant_genes)
head(annoatated_data)

# Assuming your annotated data is stored in 'annotated_data' and the second dataframe is 'other_data'

merged_data <- merge(annoatated_data, significant_genes, by.x = "ProbeName", by.y = "ProbeName")
