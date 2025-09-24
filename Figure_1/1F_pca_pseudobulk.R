library(Seurat)
library(openxlsx)
library(ggplot2)
library(scales)
library(ggbiplot)

MNP_LP <- readRDS("~/Dropbox/Overførsler/R3_CD14CD1Cprol_MNPLP.rds")
MNP_LP@images <- list()


#### F1 update to fit HM (renaming clusters ...) 22.07.04 ####
ord <- c('24','12','22','16','4','33','19','17','2','15','35','3','18','13','21','5','28','0','10','29','37','34','36','38','26','7','8','30','32','27','11','20','9','1','14','6','31','25','23')
ord <- as.character(ord)

MNP_LP@meta.data$integrated_snn_res.2.8 <- as.character(MNP_LP@meta.data$integrated_snn_res.2.8)
MNP_LP@meta.data$integrated_snn_res.2.8 <- factor(MNP_LP@meta.data$integrated_snn_res.2.8, levels=ord)
Idents(MNP_LP) <- 'integrated_snn_res.2.8'

col_pal <- hue_pal()(39)
names(col_pal) <- ord

MNP_LP@meta.data$lineage <- NA
for (cl in levels(MNP_LP@meta.data$integrated_snn_res.2.8)){
  if (cl %in% levels(MNP_LP@meta.data$integrated_snn_res.2.8)[1:3]){
    MNP_LP@meta.data[MNP_LP@meta.data$integrated_snn_res.2.8==cl,]$lineage <- "mono"}
  else if (cl %in% levels(MNP_LP@meta.data$integrated_snn_res.2.8)[4:20]){
    MNP_LP@meta.data[MNP_LP@meta.data$integrated_snn_res.2.8==cl,]$lineage <- "mac"}
  else {MNP_LP@meta.data[MNP_LP@meta.data$integrated_snn_res.2.8==cl,]$lineage <- "cDC2-like"}
}

## gene lists 
Mulder_sign <- read.xlsx("1-s2.0-S1074761321002934-mmc5.xlsx", sep = "\t", startRow = 2)

Macro <- Mulder_sign %>% filter(cluster == "Macro") %>% pull(Gene)
cMo <- Mulder_sign %>% filter(cluster == "cMO") %>% pull(Gene)
DC1 <- Mulder_sign %>% filter(cluster == "cDC1") %>% pull(Gene)
DC2 <- Mulder_sign %>% filter(cluster == "cDC2") %>% pull(Gene)
cur_lin_genes <- c(Macro, cMo, DC1, DC2)


## av pca
Idents(MNP_LP) <- MNP_LP@meta.data$integrated_snn_res.2.8
MNP_LP_av <- AverageExpression(MNP_LP, features = cur_lin_genes, assays = "RNA", return.seurat = T)
MNP_LP_av_pc <- t(MNP_LP_av$RNA$data)

leaveouts <- c()
for (cluster in rownames(MNP_LP_av_pc)){
  navs <- colnames(MNP_LP_av_pc)[which(is.na(MNP_LP_av_pc[cluster,]))]
  if (length(navs)>0){
    print(paste(cluster,navs))
    leaveouts <- append(leaveouts,navs)
  }}
saves <- colnames(MNP_LP_av_pc)[!colnames(MNP_LP_av_pc) %in% leaveouts]
MNP_LP_av_pc <- MNP_LP_av_pc[,saves]
dim(MNP_LP_av_pc) #162 genes

zero_var_cols <- apply(MNP_LP_av_pc, 2, var) == 0
MNP_LP_av_pc <- MNP_LP_av_pc[, !zero_var_cols]
MNP_pca <- prcomp(MNP_LP_av_pc, center = TRUE,scale = TRUE)
MNP_df <- as.data.frame(MNP_pca$x)
MNP_PCA_2 <- PCA(scale(MNP_LP_av_pc))

MNP_df$lineage <- NA
rownames(MNP_df) <- gsub("g", "", rownames(MNP_df))
MNP_df[rownames(MNP_df) %in% c('12','24','22'),]$lineage <- "monocytes"
MNP_df[rownames(MNP_df) %in% c('4','33','16','19'),]$lineage <- "monocytes" #"mono-ints"
MNP_df[rownames(MNP_df) %in% c('18','35','17','28','0','10','15','2','5','29','21','3','13'),]$lineage <- "macrophages"
MNP_df[rownames(MNP_df) %in% c('37','34','36'),]$lineage <- "DC2-like"#"DC-precursors"
MNP_df[rownames(MNP_df) %in% c('38'),]$lineage <- "DC2-like"#"mig-DC2"
MNP_df[rownames(MNP_df) %in% c('26','7','8','30','32','27'),]$lineage <- "DC2-like"#"DC2A"
MNP_df[rownames(MNP_df) %in% c('11','20','9','1','14','6','31','25','23'),]$lineage <- "DC2-like"#"DC2B/3"


# defining colors
col_pal_ex <- append(col_pal, c("#A39446","#5B0FEC","#0F0266"))
names(col_pal_ex) <- append(names(col_pal),c("DC2-like","macrophages","monocytes"))

rownames(MNP_pca$x) <- seq(1,39,1)
col_pal_ex_new <- col_pal_ex
names(col_pal_ex_new)[1:39] <- seq(1,39,1)
col_pal[1:20] <- "darkred"
col_pal[21:39] <- "darkblue"

ellipse_colors <- c("DC2-like" = "#A39446", "macrophages" = "#5B0FEC", "monocytes" = "#0F0266")

# plotting 
pdf("pca_pseudobulk_mulder_deg.pdf",height = 6, width = 6)
plot1 <- ggbiplot(MNP_pca, obs.scale = 1, var.scale = 1, groups = MNP_df$lineage,alpha=0, choices = c(1,2), ellipse = T, 
                  ellipse.linewidth = 0.5, ellipse.fill = T, ellipse.alpha = 0.65, circle = FALSE, var.axes = F) +
  geom_text(label = rownames(MNP_pca$x), color = col_pal, show.legend  = F, size = 3.2) +
  theme_classic()+theme(legend.title = element_text(face="bold"))
plot1 + scale_colour_manual(values = col_pal_ex_new, name = "cluster") +
  scale_fill_manual(values = ellipse_colors, name = "lineage") + 
  guides(colour = guide_legend(override.aes = list(shape=16))) + coord_fixed(ratio = 1.25) 
dev.off()












plot1 <- ggbiplot(MNP_pca,
                  obs.scale = 1, var.scale = 1,
                  groups = MNP_df$cluster,   # keep clusters for points
                  alpha = 0,
                  choices = c(1,2),
                  ellipse = FALSE,           # turn OFF built-in ellipses
                  circle = FALSE, var.axes = FALSE) +
  geom_point(aes(color = MNP_df$cluster), size = 2) +
  geom_text(label = rownames(MNP_pca$x),
            aes(color = MNP_df$cluster), show.legend = FALSE) +
  # add ellipses by lineage, with fill
  stat_ellipse(aes(x = MNP_pca$x[,1],
                   y = MNP_pca$x[,2],
                   fill = MNP_df$lineage,
                   color = MNP_df$lineage),
               geom = "polygon", alpha = 0.2, level = 0.95) +
  theme_classic() +
  theme(legend.title = element_text(face = "bold")) +
  scale_color_manual(values = col_pal, name = "Cluster") +
  scale_fill_manual(values = ellipse_colors, name = "Lineage") +
  guides(color = guide_legend(override.aes = list(shape = 16)))
plot1
