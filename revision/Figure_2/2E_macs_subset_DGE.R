##### Pseudobulk code
library(DESeq2)
library(Seurat)

# objects for pseudobulk
monomac <- readRDS("~/Downloads/monomac_traj_versionApril2022.rds")
monomac@images <- list()

monomac <- UpdateSeuratObject(monomac)
Idents(monomac) <- "tSP_clustering_F4"
macs <- subset(monomac, idents = c("M6", "M7", "M8"))


DefaultAssay(macs) <- "RNA"
seurat_obj <- macs

meta_data <- seurat_obj@meta.data
# Specify the cluster of interest (done for both M6, M7 and M8)
cluster_of_interest <- "M8"

# Create a new column to combine other clusters into a single reference group
meta_data$comparison_group <- ifelse(meta_data$tSP_clustering_F4 == cluster_of_interest, cluster_of_interest, "OtherClusters")
table(meta_data$comparison_group)
seurat_obj@meta.data <- meta_data

# Example: Create a grouping variable (e.g., combine cluster and sample)
seurat_obj$group <- paste0(seurat_obj@meta.data$patient, "_", seurat_obj@meta.data$comparison_group,"_", seurat_obj@meta.data$segment)
# Verify grouping
table(seurat_obj$group)


# Extract raw counts matrix
counts <- GetAssayData(seurat_obj, layer = "counts")
# Aggregate counts by group (cluster and sample)
pseudo_counts <- as.data.frame(counts) %>%
  t() %>%
  as.data.frame() %>%
  mutate(group = seurat_obj$group) %>%
  group_by(group) %>%
  summarise(across(everything(), sum)) %>%
  t()
# Assign group names to the pseudobulk matrix
colnames(pseudo_counts) <- pseudo_counts[1,] # First row contains group names
pseudo_counts <- pseudo_counts[-1,]          # Remove first row (now used as colnames)
# Convert to numeric
pseudo_counts <- apply(pseudo_counts, 2, as.numeric)
rownames(pseudo_counts) <- rownames(counts)
# Check the pseudobulk matrix
head(pseudo_counts)

# Create a sample info dataframe
grouping_info <- data.frame(group = colnames(pseudo_counts))
grouping_info$patient <- sapply(strsplit(as.character(grouping_info$group), "_"), `[`, 1)
grouping_info$cluster <- sapply(strsplit(as.character(grouping_info$group), "_"), `[`, 2)
grouping_info$condition <- sapply(strsplit(as.character(grouping_info$group), "_"), `[`, 3)
# Convert to factor for DESeq2
grouping_info$cluster <- factor(grouping_info$cluster)
grouping_info$patient <- factor(grouping_info$patient)
grouping_info$condition <- factor(grouping_info$condition)

# Check the metadata
head(grouping_info)


# Create a DESeq2 dataset (in 'design' can put both cluster + sample if they are not dependent of each other)
dds <- DESeqDataSetFromMatrix(countData = pseudo_counts,
                              colData = grouping_info,
                              design = ~ cluster + patient + condition)


# Run the DESeq2 differential expression analysis
dds <- DESeq(dds)


# Specify the comparison of interest in the results
# Replace "Cluster_A" and "Cluster_B" with your actual cluster names
res <- results(dds, contrast = c("cluster", cluster_of_interest, "OtherClusters"))

# Filter for significant genes based on adjusted p-value threshold
sig_genes <- as.data.frame(res) %>% filter(padj < 0.05)


# Assuming 'sig_genes' contains the filtered significant genes and 'res' has row names as gene names
sig_genes$gene <- rownames(sig_genes)  # Add gene names to the dataframe for labeling
length(sig_genes$gene)

M8_DEGs <- as.data.frame(sig_genes) %>% filter(log2FoldChange > 0)

write.xlsx(M6_DEGs, "M6_DEGs.xlsx")
write.xlsx(M7_DEGs, "M7_DEGs.xlsx")
write.xlsx(M8_DEGs, "M8_DEGs.xlsx")


# Create the Venn diagram -- Figure S2C
pdf("M6-8_VennDiagram_new.pdf")
venn <- venn.diagram(
  x = list(M6_DEGs$gene, M7_DEGs$gene, M8_DEGs$gene),
  category.names = c("M6", "M7", "M8"),
  fill = c("#59C1AE", "#54B8DD", "#4BA4F1"),
  alpha = 0.5,
  lwd = c(0, 0, 0),
  #lty = "blank",
  filename = NULL,
  output = TRUE
)
grid.draw(venn)
dev.off()

## heatmap DGEs -- Figure 2E
M6_DEGs <- M6_DEGs |> arrange(padj)
M7_DEGs <- M7_DEGs |> arrange(padj)
M8_DEGs <- M8_DEGs |> arrange(padj)
genes.to.plot <- c(M6_DEGs$gene[1:15], M7_DEGs$gene[1:15], M8_DEGs$gene[1:15])
monomactraj_Av <- AverageExpression(macs, return.seurat = T, features = genes.to.plot,
                                    group.by = c("tSP_clustering_F4", "segment"))

mycols <- rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))

pdf("2E_M6-8_top15_DGE_heatmap.pdf", height = 16)
DoHeatmap(monomactraj_Av, features = genes.to.plot, 
          draw.lines = F) + scale_fill_gradientn(colours = mycols, na.value = "white")
dev.off()



