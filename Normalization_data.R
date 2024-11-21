# Load data
raw_counts <- read.csv("E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Raw_file.csv")
gene_lengths <- read.table("Human.GRCh38.p13.annot.tsv.gz", header=TRUE, sep="\t")

gene_lengths <- read.table("Human.GRCh38.p13.annot.tsv.gz", 
                           header=TRUE, 
                           sep="\t", 
                           fill=TRUE, 
                           comment.char="#")

# Preview the GeneID column
gene_lengths <- read.csv("E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Annotation_file.csv")
head(gene_lengths$GeneID)


# Merge counts and lengths
data <- merge(raw_counts, gene_lengths[, c("GeneID", "Length")], by="row.names")
rownames(data) <- data$Row.names
data$Row.names <- NULL

 
# Calculate TPM

# Assuming 'data' has gene IDs in the first and second last columns, and gene lengths in the last column
# Extract the gene IDs (first column) and gene lengths (last column)
gene_ids <- data[, 1]  # Gene IDs in the first column
gene_lengths <- data[, ncol(data)]  # Gene lengths in the last column

# Extract raw counts (everything except the first, second last, and last column)
raw_counts <- data[, -c(1, ncol(data), ncol(data)-1)]

# Convert gene lengths from base pairs to kilobases
lengths_kb <- gene_lengths / 1000

# Calculate TPM for each gene and sample
tpm <- apply(raw_counts, 2, function(counts) {
  # Normalize by gene length (in kilobases)
  normalized_counts <- counts / lengths_kb
  # Scale to per million
  scaling_factor <- sum(normalized_counts) / 1e6
  normalized_counts / scaling_factor
})

# Convert TPM to a data frame and optionally add gene IDs
tpm_data <- as.data.frame(tpm)
rownames(tpm_data) <- gene_ids  # Assign gene IDs as row names

# View the first few rows of the TPM data
head(tpm_data)

write.csv(tpm_data, file="TPM_normalized_counts.csv", quote=FALSE)


































tpm <- apply(data[, -ncol(data)], 2, function(counts) {
  lengths_kb <- data$Length / 1000
  normalized_counts <- counts / lengths_kb
  scaling_factor <- sum(normalized_counts) / 1e6
  normalized_counts / scaling_factor
})

# Write to file
write.csv(tpm, file="TPM_normalized_counts.csv", quote=FALSE)

Raw_data <- read.csv("E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Raw_file.csv")
Metadata <- read.csv ("E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Metadata.csv")
rownames(Raw_data) <- Raw_data[,1]
Raw_data <- Raw_data[,-1]
rownames(Metadata) <- Metadata$Samples
Metadata$Samples <- NULL
all(rownames(Metadata) == colnames(Raw_data))  # Should return TRUE

# Print both rownames and column names to inspect
rownames(Metadata)
colnames(Raw_data)

# Check for mismatches
setdiff(rownames(Metadata), colnames(Raw_data))  # Samples in metadata but not in counts
setdiff(colnames(Raw_data), rownames(Metadata))  # Samples in counts but not in metadata

rownames(Metadata) <- trimws(rownames(Metadata))  # Remove leading/trailing spaces
colnames(Raw_data) <- trimws(colnames(Raw_data))  # Do the same for column names

Metadata <- Metadata[colnames(Raw_data),,drop = FALSE]
all(rownames(Metadata) == colnames(Raw_data)) 

summary(Raw_data) 
class(Raw_data)
# Check if there are any non-integer values
any(Raw_data != floor(Raw_data))  # TRUE means there are non-integer values

Metadata$Condition <- gsub(" ", "_", Metadata$Condition)
Metadata$Condition <- as.factor(Metadata$Condition)
levels(Metadata$Condition)

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = Raw_data,
                              colData = Metadata,
                              design = ~ Condition)

dds <- DESeq(dds)

res <- results(dds)
head(res)

res$padj <- p.adjust(res$pvalue, method = "BH")

significant_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
head(significant_genes)
summary(res)
write.csv(significant_genes, file = "significant_genes.csv")



cond <- unique(dds$Condition)

# Get all pairwise comparisons
pairwise_comparisons <- combn(cond, 2, simplify = FALSE)

# Store results in a list
comparison_results <- list()

for (pair in pairwise_comparisons) {
  comparison_results[[paste(pair, collapse = "_vs_")]] <- results(dds, contrast = c("Condition", pair[1], pair[2]))
}

# Print results for each comparison
comparison_results

levels(dds$Condition)
view(cond)
cond


# Create all pairwise combinations of the conditions
pairwise_comparisons <- combn(levels(dds$Condition), 2, simplify = FALSE)

# Store results in a list
comparison_results <- list()

# Loop through each pair and perform the comparison
for (pair in pairwise_comparisons) {
  # Perform the contrast between the two conditions
  comparison_results[[paste(pair, collapse = "_vs_")]] <- results(dds, contrast = c("Condition", pair[1], pair[2]))
}

# Example: View the results of astrocytomas vs oligodendrogliomas
summary(comparison_results$astrocytomas_vs_oligodendrogliomas)

# Install openxlsx if you don't have it already
install.packages("openxlsx")

# Load the package
library(openxlsx)

# Create a new workbook
wb_2 <- createWorkbook()

# Loop through each pairwise comparison and extract significant genes
for (i in seq_along(names(comparison_results))) {
  pair <- names(comparison_results)[i]
  
  # Get the results for the current comparison
  res <- comparison_results[[pair]]
  
  # Adjust p-values (if not already adjusted)
  res$padj <- p.adjust(res$pvalue, method = "BH")
  
  # Filter for significant genes (adjusted p-value < 0.05 and |log2FC| > 1)
  significant_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
  significant_genes$GeneID <- rownames(significant_genes)
  # Create a unique sheet name using a prefix and the pair name, truncated if necessary
  sheet_name <- paste0( i, "_", substr(pair, 1, 27))  # Add a unique counter
  
  # Add the significant genes to the workbook as a new sheet
  addWorksheet(wb_2, sheetName = sheet_name)
  writeData(wb_2, sheet = sheet_name, significant_genes)
}


# Save the workbook to an Excel file
saveWorkbook(wb_2, "significant_genes_comparisons.xlsx", overwrite = TRUE)
write.csv(pairwise_comparisons, "Order_of_Comparison.csv")

plotPCA(vsd, intgroup="Condition")


# Perform regularized log transformation
rlog_data <- rlog(dds, blind = FALSE)

# Perform PCA and plot
pca_plot_rlog <- plotPCA(rlog_data, intgroup = "Condition", returnData = TRUE)

# Enhance with ggplot2 for better visualization
library(ggplot2)
ggplot(pca_plot_rlog, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA of rlog Transformed Data", x = "PC1: variance explained", y = "PC2: variance explained") +
  theme(legend.position = "right")


# Use VST or rlog-transformed data
data_for_clustering <- assay(vsd)  # Replace `vsd` with your VST or rlog object

# Optionally, focus on the top variable genes for clustering
library(genefilter)
top_var_genes <- head(order(rowVars(data_for_clustering), decreasing = TRUE), 50)
data_for_clustering <- data_for_clustering[top_var_genes, ]

library(ConsensusClusterPlus)
# Input the normalized matrix (genes as rows, samples as columns)
res_cluster <- ConsensusClusterPlus(
  as.matrix(data_for_clustering),          # Input data
  maxK = 14,                                # Maximum number of clusters to test (adjust based on your data)
  reps = 1000,                             # Number of resampling iterations
  pItem = 0.8,                             # Proportion of samples to subsample in each iteration
  pFeature = 1,                            # Proportion of features (genes) to subsample
  clusterAlg = "km",                       # Clustering algorithm (e.g., hierarchical clustering)
  distance = "euclidean",                    # Distance metric (e.g., "euclidean", "pearson")
  seed = 123,                              # Set a random seed for reproducibility
  plot = "png",
  writeTable = TRUE,
  title = "Cluster top 5000"
)

# Assuming 'result' is the output of ConsensusClusterPlus
k <- 4  # Replace with the desired number of clusters (e.g., k=4)
clusters <- res_cluster[[k]]$consensusClass  # Extract the cluster assignments

cluster_assignments <- data.frame(
  SampleID = names(clusters),
  Cluster = clusters
)

Metadata_new <- read.csv ("E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Metadata.csv")
rownames(Metadata_new) <- Metadata_new[,1]
Metadata_new <- Metadata_new[colnames(Raw_data),,drop = FALSE]


metadata_with_clusters <- merge(Metadata_new, cluster_assignments, by.x = "Samples", by.y = "SampleID")

rownames(Metadata_new) <- NULL  # Optional: remove row names if no longer needed


data_for_clustering_meta <- assay(vsd)

# Transpose expression data to match sample-level metadata
transposed_expression <- t(data_for_clustering_meta)

# Convert to data frame for merging
transposed_expression <- as.data.frame(transposed_expression)
transposed_expression$SampleID <- rownames(transposed_expression)

# Merge with clinical and cluster information
expression_with_info <- merge(metadata_with_clusters, transposed_expression, by.x = "Samples", by.y = "SampleID")

# Reorder the data frame based on the Cluster column
expression_with_info <- expression_with_info[order(expression_with_info$Cluster), ]


final_expression <- expression_with_info # Remove SampleID, Condition, Cluster columns
rownames(final_expression) <- expression_with_info$Samples
final_expression <- t(final_expression)


# Create sample annotation for conditions and clusters
sample_annotation <- data.frame(
  Condition = expression_with_info$Condition,  # Replace with your condition column
  Cluster = expression_with_info$Cluster,   # Replace with your cluster assignments
  row.names = expression_with_info$Samples  # Match sample names with your columns
)

# Define colors for each condition
condition_colors <- c(
  "anaplastic astrocytomas" = "#FF6347",  # Tomato red
  "anaplastic oligodendroastrocytomas" = "#4682B4",  # SteelBlue
  "anaplastic oligodendrogliomas" = "#32CD32",  # LimeGreen
  "astrocytomas" = "#FFD700",  # Gold
  "oligodendroastrocytomas" = "#6A5ACD",  # SlateBlue
  "oligodendrogliomas" = "#8B0000",  # DarkRed
  "primary Glioblastomas" = "#FF1493",  # DeepPink
  "recurrent anaplastic astrocytomas" = "#FF8C00",  # DarkOrange
  "recurrent anaplastic oligodendroastrocytomas" = "#00CED1",  # DarkTurquoise
  "recurrent anaplastic oligodendrogliomas" = "#C71585",  # MediumVioletRed
  "recurrent astrocytomas" = "#DAA520",  # GoldenRod
  "recurrent Glioblastomas" = "#B22222",  # FireBrick
  "recurrent oligodendroastrocytomas" = "#20B2AA",  # LightSeaGreen
  "secondary Glioblastomas" = "#8A2BE2"   # BlueViolet
)

cluster_color <- c( "1" = "yellow", "2" = "green", "3" = "blue", "4" = "red")

heatmap_input <- data_for_clustering
unique(sample_annotation$Condition)
setdiff(unique(sample_annotation$Condition), names(condition_colors))
setdiff(unique(sample_annotation$Cluster), names(cluster_color))






library(ComplexHeatmap)
# Create a HeatmapAnnotation object for annotating the heatmap
ha <- HeatmapAnnotation(
  Condition = sample_annotation$Condition[sample_order],
  Cluster = sample_annotation$Cluster[sample_order],
  col = list(
    Condition = condition_colors,  # Replace with your conditions
    Cluster =  cluster_color # Customize colors for clusters
  )
)

heatmap_expression <- final_expression
heatmap_expression <- heatmap_expression[-c(1,2,3)]
heatmap_expression

heatmap <- Heatmap(
 heatmap_new_input,  # Expression data (genes x samples)
  name = "Expression",
  top_annotation = ha,  # Add condition annotations on top
  show_row_names = TRUE,  # Show gene names (you can hide it if there are too many genes)
  show_column_names = FALSE,  # Show sample names (optional)
  cluster_rows = TRUE,  # Cluster genes (rows)
  cluster_columns = FALSE,  # Cluster samples (columns)
  heatmap_legend_param = list(title = "Expression Level", legend_width = unit(4, "cm")),
  column_title = "Samples",  # Column title (for samples)
  row_title = "Genes"  # Row title (for genes)
)

draw(heatmap)
write.csv(heatmap_input_top_scaled, "E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Normalization/Heatmap_data.csv")

heatmap_input_top_scaled <- t(scale(t(heatmap_input)))
sample_order <- order(sample_annotation$Cluster)
heatmap_input_top_scaled <- heatmap_input_top_scaled[, sample_order]


heatmap_new_input <- read.csv("E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Normalization/Heatmap_data.csv")
rownames(heatmap_new_input) <- heatmap_new_input[,1]
heatmap_new_input$X <- NULL
heatmap_new_input < - as.matrix(heatmap_new_input)
class(heatmap_new_input)
str(heatmap_new_input)
heatmap_new_input <- as.matrix(heatmap_new_input)

write.csv(metadata_with_clusters, "E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Normalization/metadata_with_clusters.csv")


# Example: sample_annotation has columns Cluster and Condition
table_cluster_condition <- table(metadata_with_clusters$Cluster, metadata_with_clusters$Condition)
print(table_cluster_condition)
# Find the cluster with the most samples for each condition
max_condition_per_cluster <- apply(table_cluster_condition, 2, which.max)
print(max_condition_per_cluster)


# Extract the condition with the maximum frequency for each cluster
max_condition_per_cluster_names <- apply(table_cluster_condition, 1, function(x) names(x)[which.max(x)])

# Create a data frame for plotting
cluster_condition_df <- data.frame(
  Cluster = rownames(table_cluster_condition),
  Condition = max_condition_per_cluster_names,
  Freq = apply(table_cluster_condition, 1, max)  # Get the frequency of the most common condition in each cluster
)

# Print the resulting data frame
print(cluster_condition_df)


library(ggplot2)

library(ggplot2)

# Create the plot with a color palette and a legend
ggplot(cluster_condition_df, aes(x = Cluster, y = Freq, fill = Condition)) +
  geom_bar(stat = "identity", show.legend = TRUE) +  # Ensure the legend is shown
  labs(x = "Cluster", y = "Frequency of Most Frequent Condition", fill = "Condition") +
  theme_minimal() +
  scale_fill_manual(values = c(
    "astrocytomas" = "#FFD700", 
    "primary Glioblastomas" = "#FF1493", 
    "secondary Glioblastomas" = "#8A2BE2",
    "oligodendroastrocytomas" = "#00CED1"
    # You can add more conditions and colors as needed
  )) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability



write.csv(expression_with_info,"E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Normalization/Expression_with_info.csv")
write.csv(final_expression,"E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Normalization/final_expression.csv")

















# Load the necessary library
library(DESeq2)

# Read the CSV file without setting row.names
expression_data_raw <- read.csv("E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Cluster1_cluster2.csv", header = TRUE)

expression_data_raw$Gene_names <- rownames(expression_data_raw)
# View the first few rows to inspect the data
head(expression_data_raw)

# Check for duplicate row names
duplicated_genes <- duplicated(expression_data_raw$GeneName) | duplicated(expression_data_raw$GeneName, fromLast = TRUE)

# View rows with duplicated gene names
expression_data_raw[duplicated_genes, ]



# Read the expression data from your file (ensure you set the correct file path)
expression_data <- read.csv("E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Cluster1_cluster2.csv", header = TRUE, row.names = 1)

# Check the first few rows to make sure it's loaded correctly
head(expression_data)

# Remove rows with duplicated gene names, keeping only the first occurrence
# Read in the expression data, assuming the first column is gene names
expression_data_raw <- read.csv("E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Cluster1_cluster2.csv", header = TRUE, row.names = 1)

# If the first two rows contain metadata (e.g., sample and cluster info), remove them
expression_data_raw <- expression_data_raw[-c(1,2), ]
rownames(expression_data_raw) <- expression_data_raw[,1]

expression_data_no_duplicates <- expression_data_raw[!duplicated(expression_data_raw[, 1]), ]




# Remove rows with duplicate gene names
expression_data_no_duplicates <- expression_data_raw[!duplicated(rownames(expression_data_raw)), ]

# Verify the result
head(expression_data_no_duplicates)


write.csv(expression_data_no_duplicates,"E:/Transcriptomics/Multiomics/Transcriptomics_analysis/Normalization/expression_1_2.csv" )
# Check the result
head(expression_data_no_duplicates)


# Check if the duplicates were removed
head(expression_data_no_duplicates)






metadata_with_clusters$Cluster <- factor(metadata_with_clusters$Cluster)


dds <- DESeqDataSetFromMatrix(countData = Raw_data,
                              colData = metadata_with_clusters,
                              design = ~ Cluster)

dds <- DESeq(dds)

res <- results(dds)
head(res)


pairwise_comparisons <- combn(levels(dds$Cluster), 2, simplify = FALSE)

# Store results in a list
comparison_results <- list()

# Loop through each pair and perform the comparison
for (pair in pairwise_comparisons) {
  # Perform the contrast between the two conditions
  comparison_results[[paste(pair, collapse = "_vs_")]] <- results(dds, contrast = c("Cluster", pair[1], pair[2]))
}


wb_3 <- createWorkbook()

# Loop through each pairwise comparison and extract significant genes
for (i in seq_along(names(comparison_results))) {
  pair <- names(comparison_results)[i]
  
  # Get the results for the current comparison
  res <- comparison_results[[pair]]
  
  # Adjust p-values (if not already adjusted)
  res$padj <- p.adjust(res$pvalue, method = "BH")
  
  # Filter for significant genes (adjusted p-value < 0.05 and |log2FC| > 1)
  significant_genes <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
  significant_genes$GeneID <- rownames(significant_genes)
  # Create a unique sheet name using a prefix and the pair name, truncated if necessary
  sheet_name <- paste0( i, "_", substr(pair, 1, 27))  # Add a unique counter
  
  # Add the significant genes to the workbook as a new sheet
  addWorksheet(wb_3, sheetName = sheet_name)
  writeData(wb_3, sheet = sheet_name, significant_genes)
}


# Save the workbook to an Excel file
saveWorkbook(wb_3, "significant_genes_comparisons_new.xlsx", overwrite = TRUE)
