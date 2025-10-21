
######
#### High-dimensional profiling uncovers heterogeneity, lineage-specific precursors and inflammation-induced changes in the mononuclear phagocyte compartment of the human intestine
#### Fenton, Wulff and Väänänen et al. 
#### Code for plots in Figure 7 and S7 - Analysis of cDCs in small intestinal lamina propria of Crohn's disease patients
#### Author: Venla Väänänen 

library(Seurat)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggpubr)
library(scales)

mycols <- rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
mycols_b <- c("#bdbdbd","#d9d9d9","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")

obj6_si <- readRDS("cDC_LP_anon_VV.rds")


# Fig. 7B
DimPlot(obj6_si, group.by = 'cDC_annotations', reduction = 'umap.harmony', raster = TRUE, raster.dpi = c(702,702), pt.size = 2, cols = c("#00BFC4","#B79F00","#00BA38","#F8766D","#F564E3","#619CFF"))+ coord_fixed(ratio = 1.15)

# Fig. S7A
DimPlot(obj6_si, group.by = 'Line_labels', reduction = 'umap.harmony', alpha = 0.7, raster = TRUE, raster.dpi = c(702,702), pt.size = 2, cols = c("#00BA38","#F8766D","#B79F00","#00BFC4","#619CFF","lightgrey"))+ coord_fixed(ratio = 1.15)

# Fig. S7B
DC1score <- c("XCR1","CADM1","CLEC9A")
DC1.score <- list(DC1score)
obj6_si <- AddModuleScore(obj6_si, features = DC1.score, ctrl = 5, name = 'DC1_smallscore')

DC2score <- c("LTB","CD207","FLT3","IRF4")
DC2.score <- list(DC2score)
obj6_si <- AddModuleScore(obj6_si, features = DC2.score, ctrl = 5, name = 'DC2_smallscore')

DC3score <- c("CD1C","CD163","FCER1A","CCR7")
DC3.score <- list(DC3score)
obj6_si <- AddModuleScore(obj6_si, features = DC3.score, ctrl = 5, name = 'DC3_smallscore')

CCR7DCscore <- c('CCR7', 'LAMP3', 'CD40')
CCR7DC.score <- list(CCR7DCscore)
obj6_si <- AddModuleScore(obj6_si, features = CCR7DC.score, ctrl = 5, name = 'CCR7DC_score')

proHLAlowscore <- c("KIAA0101", "STMN1", "H2AFZ", "HIST1H4C", "HMGN2", "TUBB", "S100B", "THBS1", "CCL20")
proHLAlow.score <- list(proHLAlowscore)
obj6_si <- AddModuleScore(obj6_si, features = proHLAlow.score, ctrl = 5, name = 'proHLAlow_score')

FeaturePlot(obj6_si, features = 'DC1_smallscore1', reduction = 'umap.harmony', raster = TRUE, raster.dpi = c(702,702), pt.size = 2)+ coord_fixed(ratio = 1.15) +scale_colour_gradientn(colours=mycols)
FeaturePlot(obj6_si, features = 'DC2_smallscore1', reduction = 'umap.harmony',raster = TRUE, raster.dpi = c(702,702), pt.size = 2)+ coord_fixed(ratio = 1.15) +scale_colour_gradientn(colours=mycols)
FeaturePlot(obj6_si, features = 'DC3_smallscore1', reduction = 'umap.harmony', raster = TRUE, raster.dpi = c(702,702), pt.size = 2)+ coord_fixed(ratio = 1.15) +scale_colour_gradientn(colours=mycols)
FeaturePlot(obj6_si, features = 'CCR7DC_score1', reduction = 'umap.harmony', raster = TRUE, raster.dpi = c(702,702), pt.size = 2)+ coord_fixed(ratio = 1.15) +scale_colour_gradientn(colours=mycols)
FeaturePlot(obj6_si, features = 'proHLAlow_score1', reduction = 'umap.harmony', raster = TRUE, raster.dpi = c(702,702), pt.size = 2)+ coord_fixed(ratio = 1.15) +scale_colour_gradientn(colours=mycols)


# Fig. 7C

prop_df <- obj6_si@meta.data %>%
  group_by(Patient, Sample_Condition, cDC_annotations) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Patient, Sample_Condition) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup() %>%
  # rename patients according to map
  mutate(Patient = recode(Patient, !!!patient_map))

# Choose the cluster you want to plot
cluster_of_interest <- "proHLAlow"

plot_df <- prop_df %>%
  filter(cDC_annotations == cluster_of_interest) %>%
  select(Patient, Sample_Condition, Proportion)

long_data <- plot_df %>%
  pivot_wider(names_from = Sample_Condition, values_from = Proportion) %>%
  pivot_longer(cols = c("Non_Inflamed","Inflamed"),
               names_to = "Condition", values_to = "Proportion")

crc_df <- prop_df %>%
  filter(Patient %in% c("pat1","pat2","pat3","pat4"), 
         cDC_annotations == cluster_of_interest) %>%
  mutate(Condition = "unpaired_crc")

# Plot
ggplot() +
  # Paired patients
  geom_boxplot(data = long_data, aes(x = Condition, y = Proportion, color = Condition), width = 0.3) +
  geom_line(data = long_data, aes(x = Condition, y = Proportion, group = Patient), color = "black") +
  geom_point(data = long_data, aes(x = Condition, y = Proportion, color = Condition), size = 4) +
  # Unpaired CRC
  geom_boxplot(data = crc_df, aes(x = Condition, y = Proportion), width = 0.3, color = "black") +
  geom_point(data = crc_df, aes(x = Condition, y = Proportion), color = "black", size = 4) +
  # Stats
  stat_compare_means(data = long_data, aes(x = Condition, y = Proportion),
                     paired = TRUE, label = "p.format",
                     label.y = max(long_data$Proportion, na.rm = TRUE) + 0.02) +
  scale_y_continuous(limits = c(0, NA), labels = label_number(accuracy = 0.1)) +
  scale_color_manual(values = c("black","black")) +
  theme_classic() +
  labs(x = " ", y = "Proportion", title = cluster_of_interest)


# Fig. S7C
Idents(obj6_si) <- obj6_si@meta.data$Sample_Condition #here with crc samples
obj7 <- subset(x = obj6_si, idents = c('Inflamed','Non_Inflamed'))

meta <- obj7@meta.data
meta$cluster_condition <- paste0(meta$cDC_annotations, "_", meta$Sample_Condition)
table(meta$cluster_condition)
obj7@meta.data <- meta
obj7@meta.data

Idents(obj7) <- obj7@meta.data$cluster_condition

DefaultAssay(obj7) <- "ADT"
my_levels <- c('cDC1_Non_Inflamed','cDC1_Inflamed','cDC2_Non_Inflamed','cDC2_Inflamed','amb_Non_Inflamed','amb_Inflamed','cDC3_Non_Inflamed','cDC3_Inflamed','CCR7pos_Non_Inflamed','CCR7pos_Inflamed')
levels(obj7) <- my_levels

cluster_colors_orig <- c(
  "cDC1_Non_Inflamed" = "#F8766D",
  "cDC1_Inflamed" = "#F8766D", 
  "cDC2_Non_Inflamed" = "#B79F00",
  "cDC2_Inflamed" = "#B79F00",
  "amb_Non_Inflamed" = "#00BA38",
  "amb_Inflamed" = "#00BA38", 
  "cDC3_Non_Inflamed" = "#00BFC4",
  "cDC3_Inflamed" = "#00BFC4",
  "CCR7pos_Non_Inflamed" = "#F564E3",
  "CCR7pos_Inflamed" = "#F564E3" 
)


VlnPlot(obj7, features = 'Hu.CD1c', pt.size = 0, cols = cluster_colors_orig)
VlnPlot(obj7, features = 'Hu.CD11a', pt.size = 0, cols = cluster_colors_orig)

DefaultAssay(obj7) <- "RNA"
VlnPlot(obj7, features = 'CD207', pt.size = 0, cols = cluster_colors_orig)


# Fig. 7D

table(obj6_si$cDC_annotations)
Idents(obj6_si) <- obj6_si@meta.data$cDC_annotations
obj6_noProlif <- subset(x = obj6_si, idents = c('cDC1','cDC2','cDC3','amb','CCR7pos'))
Idents(obj6_noProlif) <- obj6_noProlif@meta.data$Sample_Condition
obj6_noProlif <- subset(x = obj6_noProlif, idents = c('Inflamed', 'Non_Inflamed'))

table(obj6_noProlif$Sample_Condition)
table(obj6_noProlif$Patient)
table(obj6_noProlif$cDC_annotations)

##### Pseudobulk code
library(DESeq2)

# objects for pseudobulk

seurat_obj <- obj6_noProlif
seurat_obj <- JoinLayers(seurat_obj) # if use seurat v5 might need this

meta_data <- seurat_obj@meta.data
# Specify the cluster of interest
cluster_of_interest <- "CCR7pos"
# Create a new column to combine other clusters into a single reference group
meta_data$comparison_group <- ifelse(meta_data$cDC_annotations == cluster_of_interest, cluster_of_interest, "OtherClusters")
table(meta_data$comparison_group)
seurat_obj@meta.data <- meta_data

# Example: Create a grouping variable (e.g., combine cluster and sample)
seurat_obj$group <- paste0(seurat_obj@meta.data$Patient, "_", seurat_obj@meta.data$comparison_group,"_", seurat_obj@meta.data$Sample_Condition)
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
grouping_info$cluster <- sapply(strsplit(as.character(grouping_info$group), "_"), `[`, 1)
grouping_info$sample <- sapply(strsplit(as.character(grouping_info$group), "_"), `[`, 2)
grouping_info$condition <- sapply(strsplit(as.character(grouping_info$group), "_"), `[`, 3)
# Convert to factor for DESeq2
grouping_info$cluster <- factor(grouping_info$cluster)
grouping_info$sample <- factor(grouping_info$sample)
grouping_info$condition <- factor(grouping_info$condition)

# Check the metadata
head(grouping_info)


# Create a DESeq2 dataset (in 'design' can put both cluster + sample if they are not dependent of each other)
dds <- DESeqDataSetFromMatrix(countData = pseudo_counts,
                              colData = grouping_info,
                              design = ~ cluster + sample + condition)


# Run the DESeq2 differential expression analysis
dds <- DESeq(dds)


# Specify the comparison of interest in the results
# Replace "Cluster_A" and "Cluster_B" with your actual cluster names
res <- results(dds, contrast = c("sample", cluster_of_interest, "OtherClusters"))

# Filter for significant genes based on adjusted p-value threshold
sig_genes <- as.data.frame(res) %>% filter(padj < 0.05)


# Assuming 'sig_genes' contains the filtered significant genes and 'res' has row names as gene names
sig_genes$gene <- rownames(sig_genes)  # Add gene names to the dataframe for labeling
length(sig_genes$gene)
# Define the gene names to label (optional: you can label only top N genes or those with high fold change)
top_genes <- sig_genes %>%
  arrange(abs(padj)) %>% # abs(padj) or desc(abs(log2FoldChange))
  head(50)  # Label top 10 genes with the highest fold change

res$gene <- rownames(res)


# Remove rows with NA in either padj or log2FoldChange to avoid plotting issues
res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]

## Set colors
keyvals <- ifelse(
  res$log2FoldChange < -0.58 & res$padj < 0.05, '#656565',
  ifelse(res$log2FoldChange > 0.58 & res$padj < 0.05, '#F564E3',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == '#656565'] <- 'OtherClusters'
names(keyvals)[keyvals == 'grey'] <- 'ns'
names(keyvals)[keyvals == '#F564E3'] <- 'CCR7pos'

# Customize color scale for conditions
color_scale <- c("CCR7pos" = "#F564E3", "OtherClusters" = "#656565", "ns" = "grey")

# Define thresholds
fc_cutoff <- 0.58
p_cutoff <- 0.05

# Plot
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = names(keyvals))) +
  geom_point(size = 2, rasterize = TRUE) +
  scale_color_manual(values = color_scale) +
  # Add horizontal and vertical threshold lines
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
  geom_text(data = top_genes, aes(label = gene), vjust = -1, hjust = 0.5, color = "black", size = 3) + #here can adjust top_genes gened or sig_genes
  # Labels and titles
  labs(
    title = "CCR7pos DEGs to other clusters",
    x = expression(Log[2]~"fold change"),
    y = expression(-Log[10]~italic(P))
  ) +
  coord_cartesian(xlim = c(-10,10), ylim = c(0,160)) + geom_point(size = 2) +
  # Customize theme for a clean look
  theme_classic() 

## show labels for selected genes

Venla_genes <- c("CCR7","LAMP3","FSCN1","CD274","CD40","PDCD1LG2","CD200","CD83","CD80","CD86")
CCR7_labels <- sig_genes %>% filter(gene %in% Venla_genes)

# Plot
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = names(keyvals))) +
  geom_point(size = 2, rasterize = TRUE) +
  scale_color_manual(values = color_scale) +
  # Add horizontal and vertical threshold lines
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
  geom_text(data = CCR7_labels, aes(label = gene), vjust = -1, hjust = 0.5, color = "black", size = 3) + 
  #here can adjust top_genes gened or sig_genes
  # Labels and titles
  labs(
    title = "CCR7pos DEGs to other clusters",
    x = expression(Log[2]~"fold change"),
    y = expression(-Log[10]~italic(P))
  ) +
  coord_cartesian(xlim = c(-10,10), ylim = c(0,160)) + geom_point(size = 2) +
  # Customize theme for a clean look
  theme_classic() 

ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = names(keyvals))) +
  geom_point(size = 2, rasterize = TRUE) +
  scale_color_manual(values = color_scale) +
  # Add horizontal and vertical threshold lines
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
  geom_text_repel(data = CCR7_labels, aes(label = gene), color = "black", size = 2.5, box.padding = 1, max.overlaps = Inf) + 
  #here can adjust top_genes gened or sig_genes
  # Labels and titles
  labs(
    title = "CCR7pos DEGs to other clusters",
    x = expression(Log[2]~"fold change"),
    y = expression(-Log[10]~italic(P))
  ) +
  coord_cartesian(xlim = c(-10,10), ylim = c(0,160)) +
  # Customize theme for a clean look
  theme_classic() 

# Fig. S7D

Idents(obj7) <- obj7@meta.data$cDC_annotations
obj_ave <- AverageExpression(obj7, assays = 'RNA', return.seurat = TRUE)

my_levels <- c('cDC1','cDC2','amb','cDC3','CCR7pos')
levels(obj_ave) <- my_levels

genes1 <- c('CD40','CD80','CD86','RELB','CD83','CD274','PDCD1LG2','CD200','FAS','ALDH1A2','SOCS1','SOCS2') # maturation and regulatory
genes2 <- c('CCR7','MYO1G','CXCL16','ADAM8','ICAM1','FSCN1','MARCKS','MARCKSL1','TLR1','TLR2','TLR3','TLR4','TLR5','TLR6','TLR7','TLR8','TLR9','MYD88','MAVS') # migration and TLRs and adaptors

DoHeatmap(obj_ave, features = c(genes1, genes2), draw.lines = FALSE)+scale_fill_gradientn(colours=mycols)
DoHeatmap(obj_ave, features = c(genes1, genes2), draw.lines = FALSE,group.colors = cluster_colors)+scale_fill_gradientn(colours=mycols)



# Fig. 7E

library(readxl)
GO_data <- read_excel(file_path)
head(GO_data)

# Add a -log10(PValue) column for better visualization
GO_data$logPValue <- -log10(GO_data$pvalue)

# Sort data by significance (-log10(PValue)) for ordered plotting
GO_data <- GO_data[order(GO_data$logPValue, decreasing = TRUE), ]

# Bar Plot: GO Terms by -log10(PValue)
ggplot(GO_data, aes(x = reorder(GOterm, logPValue), y = logPValue)) +
  geom_bar(stat = "identity", width = 0.8, fill = "#cb71ac") +  # Bar plot
  coord_flip() +  # Flip coordinates for better readability
  labs(
    title = "GO Term Enrichment",
    y = "-log10(PValue)",
    x = NULL
  ) +
  theme_test(base_size = 8)   # Minimal theme


# Fig. S7E

table(obj$cDC_annotations)
Idents(obj) <- obj@meta.data$cDC_annotations
obj6_ccr7 <- subset(x = obj, idents = c('CCR7pos'))

obj6_ccr7[["RNA"]] <- split(obj6_ccr7[["RNA"]], f = obj6_ccr7$Patient)
obj6_ccr7 <- IntegrateLayers(
  object = obj6_ccr7, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE)
ElbowPlot(obj6_ccr7, ndims = 30)
obj6_ccr7 <- RunUMAP(obj6_ccr7, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")
obj6_ccr7 <- FindNeighbors(obj6_ccr7, reduction = "harmony", dims = 1:15)
obj6_ccr7 <- FindClusters(obj6_ccr7, resolution = 0.1, cluster.name = "harmony_clusters")

DimPlot(obj6_ccr7, group.by = 'CCR7_subclusters', reduction = 'umap.harmony', raster = TRUE, raster.dpi = c(702,702), pt.size = 4, cols = c("#A23E84","#E5A6C8"))+ coord_fixed(ratio = 1.2)


# Fig. S7F

meta_data <- obj6_ccr7@meta.data

cell_counts <- meta_data %>%
  group_by(sample = orig.ident, cluster = harmony_clusters) %>%
  summarise(cell_count = n(), .groups = "drop")

cell_counts_wide <- cell_counts %>%
  pivot_wider(names_from = cluster, values_from = cell_count, values_fill = list(cell_count = 0))

cell_counts_wide <- cell_counts_wide %>%
  mutate(ratio = `2` / `1`) %>%
  # rename patients
  mutate(sample = recode(sample, !!!patient_map))

cond_df <- meta_data %>%
  distinct(orig.ident, Sample_Condition) %>%
  mutate(orig.ident = recode(orig.ident, !!!patient_map))

plot_df <- cell_counts_wide %>%
  left_join(cond_df, by = c("sample" = "orig.ident")) %>%
  select(Patient = sample, Condition = Sample_Condition, ratio)

long_data <- plot_df %>%
  pivot_wider(names_from = Condition, values_from = ratio) %>%
  pivot_longer(cols = c("Non_Inflamed","Inflamed"),
               names_to = "Condition", values_to = "Ratio")

# Plot
ggplot() +
  geom_boxplot(data = long_data, aes(x = Condition, y = Ratio, color = Condition), width = 0.3) +
  geom_line(data = long_data, aes(x = Condition, y = Ratio, group = Patient), color = "black") +
  geom_point(data = long_data, aes(x = Condition, y = Ratio, color = Condition), size = 4) +
  stat_compare_means(data = long_data, aes(x = Condition, y = Ratio),
                     paired = TRUE, label = "p.format",
                     label.y = max(long_data$Ratio, na.rm = TRUE) + 0.1) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(values = c("black","black")) +
  theme_classic() +
  labs(x = " ", y = "Cluster0/Cluster1 ratio", title = "2 vs 1 ratio")


# Fig. S7G

# Pseudobulk analysis
seurat_obj <- obj6_ccr7
seurat_obj <- JoinLayers(seurat_obj) # if use seurat v5 might need this

# Example: Create a grouping variable (e.g., combine cluster and sample)
seurat_obj$group <- paste0(seurat_obj@meta.data$Patient, "_", seurat_obj@meta.data$CCR7_subclusters)
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
grouping_info$cluster <- sapply(strsplit(as.character(grouping_info$group), "_"), `[`, 1)
grouping_info$sample <- sapply(strsplit(as.character(grouping_info$group), "_"), `[`, 2)
# Convert to factor for DESeq2
grouping_info$cluster <- factor(grouping_info$cluster)
grouping_info$sample <- factor(grouping_info$sample)

# Check the metadata
head(grouping_info)


# Create a DESeq2 dataset (in 'design' can put both cluster + sample if they are not dependent of each other)
dds <- DESeqDataSetFromMatrix(countData = pseudo_counts,
                              colData = grouping_info,
                              design = ~ cluster + sample)


# Run the DESeq2 differential expression analysis
dds <- DESeq(dds)


# Specify the comparison of interest in the results
res <- results(dds, contrast = c("sample", "2", "1"))

# Filter for significant genes based on adjusted p-value threshold
sig_genes <- as.data.frame(res) %>% filter(padj < 0.05)

#### GO TERMS FOR EACH SUBCLUSTER
# Data in a table
library(readxl)
GO_data <- read_excel(file_path)
head(GO_data)

# Add a -log10(PValue) column for better visualization
GO_data$logPValue <- -log10(GO_data$pvalue)
# Sort data by significance (-log10(PValue)) for ordered plotting
GO_data <- GO_data[order(GO_data$logPValue, decreasing = TRUE), ]

# Bar Plot: GO Terms by -log10(PValue)
ggplot(GO_data, aes(x = reorder(GOterm, logPValue), y = logPValue)) +
  geom_bar(stat = "identity", width = 0.8, fill = "#A23E84") +  # Bar plot
  coord_flip(ylim = c(0,11)) +  # Flip coordinates for better readability
  labs(
    title = "GO Term Enrichment",
    y = "-log10(PValue)",
    x = NULL
  ) +
  theme_test(base_size = 8) + theme(aspect.ratio = 1)    # Minimal theme

# Data in a table
library(readxl)
GO_data <- read_excel(file_path)
head(GO_data)

# Add a -log10(PValue) column for better visualization
GO_data$logPValue <- -log10(GO_data$pvalue)
# Sort data by significance (-log10(PValue)) for ordered plotting
GO_data <- GO_data[order(GO_data$logPValue, decreasing = TRUE), ]

# Bar Plot: GO Terms by -log10(PValue)
ggplot(GO_data, aes(x = reorder(GOterm, logPValue), y = logPValue)) +
  geom_bar(stat = "identity", width = 0.8, fill = "#E5A6C8") +  # Bar plot
  coord_flip(ylim = c(0,11)) +  # Flip coordinates for better readability
  labs(
    title = "GO Term Enrichment",
    y = "-log10(PValue)",
    x = NULL
  ) +
  theme_test(base_size = 8) + theme(aspect.ratio = 0.5) # Minimal theme



# Fig. 7G

## tSpace script from another file
# load in the object and add CCR7+ DC subset anotations

DimPlot(obj, group.by = 'cDC_annotations', reduction = 'umap3D.tpc', raster = TRUE, raster.dpi = c(702,702), pt.size = 2, cols = c("#00BFC4","#B79F00","#00BA38","#F8766D","#cb71ac"))+ coord_fixed(ratio = 1)
DimPlot(obj8, group.by = 'CCR7_Subclusters', reduction = 'umap3D.tpc', raster = TRUE, raster.dpi = c(702,702), pt.size = 3, cols = c("#A23E84","#E5A6C8"))+ coord_fixed(ratio = 1)


# Fig. 7H

# split object only to CD samples and then by cluster

table(obj6_si$Sample_Condition)
Idents(obj6_si) <- obj6_si@meta.data$Sample_Condition
obj6_CD <- subset(x = obj6_si, idents = c('Inflamed', 'Non_Inflamed'))

table(obj6_CD$cDC_annotations)
Idents(obj6_CD) <- obj6_CD@meta.data$cDC_annotations
obj6_cDC1 <- subset(x = obj6_CD, idents = c('cDC1'))
obj6_cDC2 <- subset(x = obj6_CD, idents = c('cDC2'))
obj6_amb <- subset(x = obj6_CD, idents = c('amb'))
obj6_cDC3 <- subset(x = obj6_CD, idents = c('cDC3'))
obj6_CCR7pos <- subset(x = obj6_CD, idents = c('CCR7pos'))


##### Pseudobulk code
library(DESeq2)

# objects for pseudobulk

seurat_obj <- obj6_CCR7pos
seurat_obj <- JoinLayers(seurat_obj) # if use seurat v5 might need this

# Example: Create a grouping variable (e.g., combine cluster and sample)
seurat_obj$group <- paste0(seurat_obj@meta.data$Patient, "_", seurat_obj@meta.data$Sample_Condition)
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
grouping_info$cluster <- sapply(strsplit(as.character(grouping_info$group), "_"), `[`, 1)
grouping_info$sample <- sapply(strsplit(as.character(grouping_info$group), "_"), `[`, 2)
# Convert to factor for DESeq2
grouping_info$cluster <- factor(grouping_info$cluster)
grouping_info$sample <- factor(grouping_info$sample)

# Check the metadata
head(grouping_info)


# Create a DESeq2 dataset (in 'design' can put both cluster + sample if they are not dependent of each other)
dds <- DESeqDataSetFromMatrix(countData = pseudo_counts,
                              colData = grouping_info,
                              design = ~ cluster + sample)


# Run the DESeq2 differential expression analysis
dds <- DESeq(dds)


# Specify the comparison of interest in the results
res <- results(dds, contrast = c("sample", "Inflamed", "Non"))

# Filter for significant genes based on adjusted p-value threshold
sig_genes <- as.data.frame(res) %>% filter(padj < 0.05)


# Assuming 'sig_genes' contains the filtered significant genes and 'res' has row names as gene names
sig_genes$gene <- rownames(sig_genes)  # Add gene names to the dataframe for labeling
length(sig_genes$gene)
# Define the gene names to label (optional: you can label only top N genes or those with high fold change)
top_genes <- sig_genes %>%
  arrange(abs(padj)) %>% # abs(padj) or desc(abs(log2FoldChange))
  head(30)  # Label top 10 genes with the highest fold change

res$gene <- rownames(res)


# Remove rows with NA in either padj or log2FoldChange to avoid plotting issues
res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]

## Set colors
keyvals <- ifelse(
  res$log2FoldChange < -0.58 & res$padj < 0.05, '#1C76BC',
  ifelse(res$log2FoldChange > 0.58 & res$padj < 0.05, '#EC2227',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == '#1C76BC'] <- 'Non'
names(keyvals)[keyvals == 'grey'] <- 'ns'
names(keyvals)[keyvals == '#EC2227'] <- 'Inflamed'

# Customize color scale for conditions
color_scale <- c("Inflamed" = "#EC2227", "Non" = "#1C76BC", "ns" = "grey")

# Define thresholds
fc_cutoff <- 0.58
p_cutoff <- 0.05

# Plot
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = names(keyvals))) +
  geom_point(size = 2) +
  scale_color_manual(values = color_scale) +
  # Add horizontal and vertical threshold lines
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
  #geom_text(data = top_genes, aes(label = gene), vjust = -1, hjust = 0.5, color = "black", size = 3) + #here can adjust top_genes gened or sig_genes
  # Labels and titles
  labs(
    title = "CCR7pos",
    x = expression(Log[2]~"fold change"),
    y = expression(-Log[10]~italic(P))
  ) +
  coord_cartesian(xlim = c(-5,5), ylim = c(0,5)) + 
  # Customize theme for a clean look
  theme_classic() 




# Fig. 7I

IFNscore <- c("STAT1",
              "IRF7",
              "CXCL9",
              "CXCL10",
              "CXCL11",
              "TNFSF10",
              "IFITM1",
              "IFITM2",
              "IFITM3")
IFN.score <- list(IFNscore)
obj8 <- AddModuleScore(obj8, features = IFN.score, ctrl = 5, name = 'IFN_score')

FeaturePlot(obj8, features = 'IFN_score1', reduction = "umap3D.tpc", raster = TRUE, raster.dpi = c(702,702), pt.size = 3) +scale_colour_gradientn(colours=mycols_b)+ coord_fixed(ratio = 1)














































