### Creating seurat object of APCs ###
library("Seurat")
library("SeuratObject")
library("SeuratDisk")
library("SeuratWrappers")
library("UCell")
library("openxlsx")
source("/home/projects/dtu_00062/people/camlem/atlas/integration_helper_fnc.R")

# save date for file naming
date <- format(Sys.time(), "%Y_%m_%d")

# file paths
file_path <- "/home/projects/dtu_00062/people/camlem/data/02_filtering/"
out <- "/home/projects/dtu_00062/people/camlem/APC_scrnaseq/results"
dir.create(out)

# loading metadata
metadata <- readRDS("/home/projects/dtu_00062/people/camlem/APC_scrnaseq/metadata.rds")

# extracting sample names
samples <- metadata$sample_ID

# making a named list of rds objects
files <- list.files(file_path)
names <- gsub("\\.rds|\\.RDS|\\.Rdata|\\.rdata", "", files)
files <- paste0(file_path, files)
names(files) <- names
files <- files[samples]

# reading data
integrated_data <- map(files, readRDS)

# rename cells using object names as prefix to avoid duplicate names for cells
for (i in names(integrated_data)) {
  integrated_data[[i]] <- RenameCells(integrated_data[[i]],
                                      add.cell.id = i
  )
}

# adding metadata
integrated_data <- lapply(integrated_data, add_meta_data, metadata,
                          out, date)

# merging data 
merged_obj <- merge(
  x = integrated_data[[1]],
  y = integrated_data[2:length(integrated_data)],
)

# removing unused assays decontX + SCT
DefaultAssay(merged_obj) <- "RNA"
merged_obj[["decontX"]] <- NULL
merged_obj[["SCT"]] <- NULL

# joining layers
merged_obj <- JoinLayers(merged_obj)

# splitting by sequencing batch
merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$Sequencing_Batch)

# processing
merged_obj <- merged_obj %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

# clustering
merged_obj <- FindNeighbors(merged_obj, dims = 1:30, reduction = "pca")
merged_obj <- FindClusters(merged_obj, resolution = 2, cluster.name = "unintegrated_clusters_res.2")

merged_obj <- RunUMAP(merged_obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# visualize by batch and clusters
plot <- DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = c("Sequencing_Batch", "Patient"))
ggsave(filename = paste0(out, "/umaps_unintegrated.png", sep = ""), plot = plot, 
       width = 12, height = 8, dpi = 300, units = "in")


# integration
merged_obj <- IntegrateLayers(
  object = merged_obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

# clustering 
merged_obj <- FindNeighbors(merged_obj, reduction = "harmony", dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 2, cluster.name = "harmony_clusters")

merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

plot <- DimPlot(
  merged_obj,
  reduction = "umap.harmony",
  group.by = c("Sequencing_Batch", "Patient", "Disease", "harmony_clusters"),
  ncol = 2
)
ggsave(filename = paste0(out, "/umaps_harmony_integrated.png", sep = ""), plot = plot, 
       width = 12, height = 8, dpi = 300, units = "in")



# re-joining layers
merged_obj <- JoinLayers(merged_obj)

# plotting expression
hla.score <- list(c("HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DOB", 
                    "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "HLA-DRB1"))
pDC.genes <- c("IL3RA", "CLEC4C", "GZMB")
other <- c("CD3E", "CD79A", "VWF", "MS4A2", "COL3A1", "NRXN1")
genes.to.plot <- c("ITGAX", "CLEC9A", "CD1C", "CD14", "FCGR3A", "LAMP3")

merged_obj <- AddModuleScore(merged_obj, features = hla.score, name = "HLA.score")
merged_obj <- AddModuleScore_UCell(merged_obj, features = hla.score, name = "HLA.score")

p1 <- FeaturePlot(merged_obj, reduction = "umap.harmony", features = pDC.genes, ncol = 3)
p2 <- FeaturePlot(merged_obj, reduction = "umap.harmony", features = other, ncol = 3)
p3 <- FeaturePlot(merged_obj, reduction = "umap.harmony", features = "HLA.score1") + scale_colour_viridis(option = "magma")
p4 <- FeaturePlot(merged_obj, reduction = "umap.harmony", features = genes.to.plot, ncol = 3)

ggsave(filename = paste0(out, "/pDCs_harmony_integrated.png", sep = ""), plot = p1, 
       width = 12, height = 8, dpi = 300, units = "in")

ggsave(filename = paste0(out, "/other_cells_harmony_integrated.png", sep = ""), plot = p2, 
       width = 12, height = 8, dpi = 300, units = "in")

ggsave(filename = paste0(out, "/HLA_score_harmony_integrated.png", sep = ""), plot = p3, 
       width = 12, height = 8, dpi = 300, units = "in")

ggsave(filename = paste0(out, "/gene_plots_harmony_integrated.png", sep = ""), plot = p4, 
       width = 12, height = 8, dpi = 300, units = "in")


# differential gene expression for cluster markers
Idents(merged_obj) <- "seurat_clusters"
markers <- FindAllMarkers(merged_obj, only.pos = T, logfc.threshold = 0.25)
markers <- markers %>%
  group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.xlsx(markers, "/home/projects/dtu_00062/people/camlem/APC_scrnaseq/cluster_markers_top50.xlsx")

saveRDS(merged_obj, paste0("/home/projects/dtu_00062/people/camlem/APC_scrnaseq/", 
                           date, "harmony_integrated_obj.rds", sep = ""))


### subsetting to HLA high populations ###
Idents(merged_obj) <- "seurat_clusters"
plot <- DimPlot(merged_obj, reduction = "umap.harmony") + NoLegend()
LabelClusters(plot = plot, id = "ident")

clusters.keep <- c("28", "26", "20", "13", "48", "2", "46", "9", "31", "25", "8", "11", "12", "30", "5", "1", "22",
                  "0", "16", "15", "4", "33", "10", "3", "19", "45", "14", "41", "39", "38")

sub <- subset(merged_obj, idents = clusters.keep)
DimPlot(sub, reduction = "umap.harmony") + NoLegend()

# re-process
sub[["RNA"]] <- split(sub[["RNA"]], f = sub$Sequencing_Batch)
sub <- sub %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

# integration
sub <- IntegrateLayers(
  object = sub, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

# clustering 
sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:30)
sub <- FindClusters(sub, resolution = 2, cluster.name = "harmony_clusters")

sub <- RunUMAP(sub, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

# plotting
plot <- DimPlot(
  sub,
  reduction = "umap.harmony",
  group.by = c("Sequencing_Batch", "Patient", "Disease", "harmony_clusters"),
  ncol = 2
)
ggsave(filename = paste0(out, "/sub_umaps_harmony_integrated.png", sep = ""), plot = plot, 
       width = 12, height = 8, dpi = 300, units = "in")


plot <- DimPlot(
  sub,
  reduction = "umap.harmony",
  group.by = c("Phase", "harmony_clusters"),
  ncol = 2
)
ggsave(filename = paste0(out, "/sub_cc_umaps_harmony_integrated.png", sep = ""), plot = plot, 
       width = 12, height = 8, dpi = 300, units = "in")


# re-joining layers
sub <- JoinLayers(sub)

# plotting expression
hla.score <- list(c("HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DOB", 
                    "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "HLA-DRB1"))
pDC.genes <- c("IL3RA", "CLEC4C", "GZMB")
other <- c("CD3E", "CD79A", "VWF", "MS4A2", "COL3A1", "NRXN1")
genes.to.plot <- c("ITGAX", "CLEC9A", "CD1C", "CD14", "FCGR3A", "LAMP3")

sub <- AddModuleScore(sub, features = hla.score, name = "HLA.score.sub")
#sub <- AddModuleScore_UCell(sub, features = hla.score, name = "HLA.score")

p1 <- FeaturePlot(sub, reduction = "umap.harmony", features = pDC.genes, ncol = 3)
p2 <- FeaturePlot(sub, reduction = "umap.harmony", features = other, ncol = 3)
p3 <- FeaturePlot(sub, reduction = "umap.harmony", features = "HLA.score.sub1") + scale_colour_viridis(option = "magma")
p4 <- FeaturePlot(sub, reduction = "umap.harmony", features = genes.to.plot, ncol = 3)

ggsave(filename = paste0(out, "/sub_pDCs_harmony_integrated.png", sep = ""), plot = p1, 
       width = 12, height = 8, dpi = 300, units = "in")

ggsave(filename = paste0(out, "/sub_other_cells_harmony_integrated.png", sep = ""), plot = p2, 
       width = 12, height = 8, dpi = 300, units = "in")

ggsave(filename = paste0(out, "/sub_HLA_score_harmony_integrated.png", sep = ""), plot = p3, 
       width = 12, height = 8, dpi = 300, units = "in")

ggsave(filename = paste0(out, "/sub_gene_plots_harmony_integrated.png", sep = ""), plot = p4, 
       width = 12, height = 8, dpi = 300, units = "in")


# differential gene expression for cluster markers
Idents(sub) <- "harmony_clusters"
markers <- FindAllMarkers(sub, only.pos = T, logfc.threshold = 0.25)
markers <- markers %>%
  group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

write.xlsx(markers, "/home/projects/dtu_00062/people/camlem/APC_scrnaseq/sub_cluster_markers_top50.xlsx")

saveRDS(sub, paste0("/home/projects/dtu_00062/people/camlem/APC_scrnaseq/", 
                           date, "_sub_harmony_integrated_obj.rds", sep = ""))



