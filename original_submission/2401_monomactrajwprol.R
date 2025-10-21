#########################################################
######### V1 - PCA as input for tspace ##################
#########################################################

setwd("/Volumes/Mucosal-Immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/trajectories/tspace/monomac_wprol")
MNP_LP <- readRDS("/Volumes/Mucosal-Immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R3_CD14CD1Cprol_MNP_LP.rds")

## Read in MNP LP R3
########### After supercomputer part #####
#### data read in and variables ####
### Clear and read in data ###
rm(list=ls())

### variables
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
# colour string for imputation and overlays
mycols_b <- c("#bdbdbd","#d9d9d9","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
# project
project <- "MNP_LP"
#trajectory subset
traj_sub <- "monomac_woprol"

##### Curated_lineage genes #####
cur_lin_genes <- c('S100A9','S100A8','FCN1','VCAN','SERPINA1','APOBEC3A','SLC11A1','S100A6','S100A4','FPR1','PLCB1','AQP9','SLC2A6','S100A12','DNAAF1','THBS1','TNFAIP6','CLEC4E','SOD2','CXCL8','IL1B','PSAP','FTL','MS4A7','IER3','FOLR2','CTSB','KLF2','SEPP1','C1QA','C1QB','SLC40A1','MAF','DAB2','SLCO2B1','PDK4','C1QC','STAB1','GADD45G','CTSD','C1QB','SEPP1','FUCA1','APOE','C1QC','C1QA','TMEM176B','LGALS3','PLA2G7','PRDX1','CTSC','TMEM176A','PLD3','SERPINF1','DUSP1','CD207','CST7','DSG2','PKIB','PPP1R14A','PPA1','CCL22','SLC38A1','LTB','DUSP5','FAM110A','HIC1','RALA','GNG7','HLA-DQB1','PLAC8','TRAF5','NAP1L1','G3BP2','CD52','ENHO','AIM1','ICAM3','CREM','SPIB','TUBA1A','NRARP','TCTN3','CST7','ARL4C','EZR','UPF2','ADAM8','IL1R2','CD1C','VDR','IFITM2','TRMT6','SDPR','C1orf162','XCR1','CLEC9A','S100B','SOX4','CCR7','LAMP3','CD40', 'ZBTB16','HEY1','MYC','KLF7','NFKBIA','NFKBIZ','MAFB','BHLHE41','NR1H3','MAF','CEBPA','ID3','NR2F2','BATF3','GTF2B','HIC1','IRF4','ZEB1','ZNF165', 'FLT3')

### data from supercomputer
subs1 <- readRDS(file="/Volumes/LWulffExD/Projects/FentonWulffData/R6/trajectories/tspace/monomac_wprol/MNP_LP_tspacefile.rds")
# head(subs1$ts_file)
# head(subs1$umap_embbeding$layout)
# head(subs1$umap_embbeding$data)
# head(subs1$umap_embbeding$knn$indexes)
# head(subs1$umap_embbeding$knn$distances)
# head(subs1$umap_embbeding$config)
# head(subs1$tspace_matrix)

#### vsulaization df ####
visu <- subs1$ts_file
monomac_subs <- subset(MNP, cells = rownames(visu))
visu <- cbind(visu, monomac_subs@meta.data[rownames(visu),])
visu <- visu[!visu$F2F5 %in% c("cDC2","cDC3","amb","HLA low"),]
#visu <- cbind(visu, DC_subs@reductions$umap@cell.embeddings[rownames(visu),])
#visu <- cbind(visu, idents=Idents(DC_subs)[rownames(visu)])


#distance from monocytes
#mono.trajectories <- subs1$tspace_matrix[,which(colnames(subs1$tspace_matrix) %in% paste0('T_', visu[which(visu$integrated_snn_res.2.8=='24'), 'Index']))]
#colnames(mono.trajectories) <- c("T_1","T_2","T_3","T_4","T_5")
#visu <- cbind(visu, dist_mono=mono.trajectories)
#visu$T_mean <-rowMeans(mono.trajectories)

#### Adding new trajectory based clustering ####
visu_obj <- subset(monomac_subs, cells = rownames(visu))
visu_obj@reductions$pca@cell.embeddings <- as.matrix(visu[,startsWith(colnames(visu),"PC")])

ElbowPlot(visu_obj)


#visu_obj <- FindNeighbors(visu_obj, reduction = "pca", dims = 1:15)
#res <- seq(0, 1, 0.1)
#visu_obj <- FindClusters(visu_obj, resolution = res)
#res <- seq(1, 2.8, 0.3)
#visu_obj <- FindClusters(visu_obj, resolution = res)
#clustree(visu_obj)

# new_clusters <- visu_obj@meta.data[,startsWith(colnames(visu_obj@meta.data),"int")]
# colnames(new_clusters) <- str_replace(colnames(new_clusters),"integrated_snn","traj")
# visu <- cbind(visu, new_clusters)

visu_obj <- RunUMAP(visu_obj, dims = 1:15, n.components = 2)
DimPlot(visu_obj, group.by = "F2F5", label = T)
visu <- visu[,!startsWith(colnames(visu),"UMAP_")]
visu <- cbind(visu, visu_obj@reductions$umap@cell.embeddings[rownames(visu),])

visu_obj@meta.data$F2F5 <- factor(visu_obj@meta.data$F2F5, levels = c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M11","mono/mac"))

png(paste(dato,traj_sub,"res.2.8.png",sep = "_"),height = 750, width = 750, res =150)
ggplot(visu, aes(x=UMAP_1,y=UMAP_2))+geom_point(size=0.3,colour="lightgrey")+
  geom_point(data = visu[visu$F2F5=="M11",],
             aes(x=UMAP_1,y=UMAP_2), colour="green")+
  geom_point(data = visu[visu$F2F5=="mono/mac",],
             aes(x=UMAP_1,y=UMAP_2), colour="red")+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()

VlnPlot(visu_obj, group.by = "F2F5", features = c("ADAMDEC1","DNASE1L3","CD63"),pt.size = 0, ncol = 2)+NoLegend()

png(paste(dato,traj_sub,"res.2.8_3D_UMAP1&2.png",sep = "_"),height = 750, width = 750, res =150)
ggplot(visu, aes(x=UMAP_1,y=UMAP_2,color=res.2.8))+geom_point(size=0.3)+theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()
png(paste(dato,traj_sub,"res.2.8_3D_UMAP1&3.png",sep = "_"),height = 750, width = 750, res =150)
ggplot(visu, aes(x=UMAP_1,y=UMAP_3,color=res.2.8))+geom_point(size=0.3)+theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()
png(paste(dato,traj_sub,"res.2.8_3D_UMAP2&3.png",sep = "_"),height = 750, width = 750, res =150)
ggplot(visu, aes(x=UMAP_2,y=UMAP_3,color=res.2.8))+geom_point(size=0.3)+theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()

ggplot(visu[!is.na(visu$branch),], aes(x=T_mean*3, y=branch,colour=res.2.8))+geom_point(position = "jitter", size=0.3)+theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))

ggplot(visu[!is.na(visu$branch),], aes(x=T_mean*3, y=branch,colour=res.2.8))+geom_point(position = "jitter", size=0.3)+theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_viridis(option=magma)


png(paste(dato,traj_sub,"int_res.2.5_Tomscolouring.png",sep = "_"),height = 750, width = 750, res =150)
ggplot(visu, aes(x=umap1,y=umap2,color=idents))+geom_point(size=0.3)+theme_classic()+scale_color_manual(values = col_2.5)+
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()
png(paste(dato,traj_sub,"ccPhaseg.png",sep = "_"),height = 750, width = 750, res =150)
ggplot(visu, aes(x=umap1,y=umap2,color=Phase))+theme_classic()+geom_point(size=0.3)+theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()
png(paste(dato,traj_sub,"Patient.png",sep = "_"),height = 750, width = 750, res =150)
ggplot(visu, aes(x=umap1,y=umap2,color=patient))+theme_classic()+geom_point(size=0.3)+theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()
png(paste(dato,traj_sub,"Pseudot_mono.png",sep = "_"),height = 750, width = 750, res =150)
ggplot(visu, aes(x=tPC1,y=tPC4,color=T_mean))+
  geom_point(size=0.3)+theme_classic()+theme_classic()+
  scale_color_gradientn(colours = c('magenta', 'gold', 'black'))+
  labs(color="pseudotime")+ggtitle("Pseudotime from monocytes")
dev.off()
for (res in colnames(new_clusters)){
  plot <- ggplot(visu, aes(x=umap1,y=umap2,color=visu[,res]))+theme_classic()+geom_point(size=0.5)+theme_classic()+
    guides(colour = guide_legend(override.aes = list(size=3)))+labs(color = res)
  
  png(paste(dato,traj_sub,res,".png",sep = "_"),height = 750, width = 750, res =150)
  print(plot)
  dev.off()}

#### 3Ds ####
visu$col0.5 <- NA
for (i in 1:length(unique(visu$integrated_snn_res.0.5))){
  cl <- sort(unique(visu$integrated_snn_res.0.5))[i]
  visu[visu$integrated_snn_res.0.5 == cl,]$col0.5 <- hue_pal()(length(unique(visu$integrated_snn_res.0.5)))[i]
  print(hue_pal()(length(unique(visu$integrated_snn_res.0.5)))[i])
  print(cl)
}

plot3d(x=visu[,"umap1"],y=visu[,"umap2"],z=visu[,"tPC4"],
       col=visu$col0.5,
       xlab = "UMAP_1",ylab = "UMAP_2",zlab = "tPC4")
writeWebGL(dir ="/Volumes/Mucosal-Immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/trajectories/tspace/monomac_woprol/",filename=paste(dato,project,traj_sub,"res2.8_3D.html", sep="_"))


#### Plots with gene exp ####
#KIAA0101, FCGR3A, and CD1C/CD207?
#genes <- c('CD14','CD55','CSF1R','CD1C','S100A9','FCN1','CD9','SOD2','MRC1','CD163','ITGAX','ITGAL','LYVE1','LYZ','CD209','MMP9','MMP12','CXCL8','CCL3','JUN','MAF','KIA0101')
genes <- cur_lin_genes #c("CLEC9A","XCR1","CD1C","CD207","CLEC12A","S100B","SOX4","IL3RA","LILRB4","IL25RA","CCR7","CD40L","LAMP3","KIAA0101")
for (gene in genes){
  if (gene %in% rownames(DC_subs@assays$integrated@data)){
    gene_exp <- DC_subs@assays$integrated@data[gene,rownames(visu[!is.na(visu$branch),])]}
  else if (gene %in% rownames(DC_subs@assays$RNA@data)){
    gene_exp <- DC_subs@assays$RNA@data[gene,rownames(visu[!is.na(visu$branch),])]}
  if (gene %in% rownames(DC_subs@assays$RNA@data)){
    plot <- ggplot(visu[!is.na(visu$branch),], aes(x=T_mean*3, y=branch,colour=gene_exp))+geom_point(position = "jitter", size=0.3)+theme_classic()+
      scale_color_viridis_c(option="inferno", name=paste(gene))+xlab("pseudotime")
    png(paste(dato,traj_sub,"overlay",gene,".png",sep = "_"),height = 750, width = 750, res =150)
    print(plot)
    dev.off()
    # plot2 <- ggplot(visu, aes(x=Mono.trajectories[,3],y=gene_exp,color=idents))+geom_point(size=0.3)+
    #   theme_classic()+scale_color_manual(values = col_2.5)+labs(x="Peudotime from monocytes")+
    #   guides(colour = guide_legend(override.aes = list(size=3)))
    # png(paste(dato,traj_sub,"pseudotime",gene,"Tomcolouring.png",sep = "_"),height = 750, width = 750, res =150)
    # print(plot2)
    # dev.off()
  }
  else {print(paste(gene,"not in data set."))}
}

genes <- c("CLEC9A","XCR1","CD1C","CD207","CLEC12A","S100B","SOX4","IL3RA","LILRB4","IL25RA","CCR7","CD40L","LAMP3","KIAA0101","CD5")
for (gene in genes){
  if (gene %in% rownames(DC_subs@assays$integrated@data)){
    gene_exp <- DC_subs@assays$integrated@data[gene,rownames(visu[!is.na(visu$branch) & visu$orig.ident=="cLP_pat4",])]}
  else if (gene %in% rownames(DC_subs@assays$RNA@data)){
    gene_exp <- DC_subs@assays$RNA@data[gene,rownames(visu[!is.na(visu$branch) & visu$orig.ident=="cLP_pat4",])]}
  if (gene %in% rownames(DC_subs@assays$RNA@data)){
    plot <- ggplot(visu[!is.na(visu$branch) & visu$orig.ident=="cLP_pat4",], aes(x=T_mean*3, y=branch,colour=gene_exp))+geom_point(position = "jitter", size=0.3)+theme_classic()+
      scale_color_viridis(option="inferno")
    png(paste(dato,traj_sub,"overlay",gene,".png",sep = "_"),height = 750, width = 750, res =150)
    print(plot)
    dev.off()
    # plot2 <- ggplot(visu, aes(x=Mono.trajectories[,3],y=gene_exp,color=idents))+geom_point(size=0.3)+
    #   theme_classic()+scale_color_manual(values = col_2.5)+labs(x="Peudotime from monocytes")+
    #   guides(colour = guide_legend(override.aes = list(size=3)))
    # png(paste(dato,traj_sub,"pseudotime",gene,"Tomcolouring.png",sep = "_"),height = 750, width = 750, res =150)
    # print(plot2)
    # dev.off()
  }
  else {print(paste(gene,"not in data set."))}
}

#### Plotting in heatmap ####
visu_obj@meta.data$idents <- visu$idents
Idents(visu_obj) <- visu_obj@meta.data$idents
visu_obj@meta.data$pseudotime <- Mono.trajectories[,3]

DoMultiBarHeatmap(visu_obj, features = genes,
                  group.by = "integrated_snn_res.0.4", additional.group.by = c("idents"))+NoLegend()
#additional.group.sort.by = c("integrated_snn_res.0.2"),cols.use = )#+ NoLegend() 



##### Can I do a better umap? ####
library(umap) #R 3.5
ggplot(visu, aes(x=umap1,y=umap2,color=integrated_snn_res.2.8))+geom_point(size=0.3)+theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = col_2.8)

umap.conf <- umap.defaults
umap.conf$n_neighbors <- 10
umap.conf$metric <- 'euclidean'
umap.conf$min_dist <- 0.01
umap.conf$n_components <- 3
set.seed(1111)
umap.ts <- umap::umap(subs1$ts_file[,startsWith(colnames(subs1$ts_file), "PC")], config = umap.conf)
colnames(umap.ts$layout) <- c('New_umap1', 'New_umap2',"New_umap3")
visualization <- cbind(visu, umap.ts$layout)
colnames(visualization)

ggplot(visualization, aes(x=New_umap1,y=New_umap2,color=integrated_snn_res.2.8))+geom_point(size=0.3)+theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))#+scale_color_manual(values = col_2.8)



### save as csv for JMP usage ####
write.csv(visualization,file="mono_traj_v1.csv")

visu_obj <- subset(monomac_subs, cells = rownames(visu))
visu_obj@reductions$pca@cell.embeddings <- as.matrix(visu[,startsWith(colnames(visu),"PC")])
colnames(umap.ts$layout) <- c('UMAP_1', 'UMAP_2',"UMAP_3")
visu_obj@reductions$umap@cell.embeddings <- as.matrix(umap.ts$layout)
DimPlot(visu_obj, group.by="integrated_snn_res.2.8", dims=c(1,3))+scale_color_manual(values = col_2.8)
visu_obj@meta.data <- visu_obj@meta.data[,c(1:8,21:22,42)]
colnames(visu_obj@meta.data)[11] <- "res.2.8"
visu_obj@meta.data <- cbind(visu_obj@meta.data, pseudotime=visualization$T_mean)
ElbowPlot(visu_obj)

visu_obj <- FindNeighbors(visu_obj, reduction = "pca", dims = 1:15)
res <- seq(0, 1, 0.1)
visu_obj <- FindClusters(visu_obj, resolution = res)

head(visu_obj@meta.data)

FeaturePlot(visu_obj, features = "pseudotime", dims=c(1,3))
saveRDS(visu_obj, file="/Volumes/Mucosal-Immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R3_CD14CD1C_monomac_traj_noprol.rds")

visu_obj <- readRDS("/Volumes/Mucosal-Immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R3_CD14CD1C_monomac_traj_noprol.rds")

colnames(visu_obj@meta.data)
DimPlot(visu_obj, group.by="res.2.8", c(1,3))
