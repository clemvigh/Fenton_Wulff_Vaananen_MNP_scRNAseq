###################################################################################
###### ----------------- normalising CITE-seq with DSB ----------------- ##########
###################################################################################
setwd("/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/citeseq")

rm(list=ls())

#### Read in libraries ####
library(gplots)
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
library(ggplot2)
library(viridis)
library(dsb)
library(Seurat)

#### Variables to use ####
#date in format year_month_day
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
# project
project <- "citeseq_MNP"

#### Read in raw data to identify empty droplets/debris ####
## EXAMPLE OF DSB NORMALIZATION FOR ONE SAMPLE 
## ALL OF THESE ARE CONCATENATED AND THEN ADDED TO SEURAT OBJECTS
## (SEE LINE 96)
pat04_SILP_GEX <- Read10X(data.dir = "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/17/filtered_feature_bc_matrix/")
pat04_SILP_CITE <- Read10X(data.dir = "/Volumes/Mucosal-immunology/WA group/Tom and Line data/patient_4/17/citeseq_data/umi_count/",gene.column=1)
rownames(pat04_SILP_CITE) <- c("Iso_IgG1","Iso_IgG2a","Iso_IgG2b","CD14","CD209","CD1c","CD55","CD5","CD207","CD206","CD11a","CD103","CD11c","unmapped")
colnames(pat04_SILP_CITE) <- paste0(colnames(pat04_SILP_CITE),"-1")
dim(pat04_SILP_GEX)
dim(pat04_SILP_CITE)

###### As Seurat objects ######

#### pat4_SILP ####
pat4_SILP <- CreateSeuratObject(counts = pat04_SILP_GEX, project = project, min.cells = 3, min.features = 100)
pat4_SILP[["percent.mt"]] <- PercentageFeatureSet(pat4_SILP, pattern = "^MT-")

GEX_cells <- rownames(pat4_SILP@meta.data)
CITE_cells <- colnames(pat04_SILP_CITE)
dim(pat04_SILP_GEX)
#dim(pat04_SILP_CITE)
length(GEX_cells[GEX_cells %in% CITE_cells]) #how many cells are we excluding from the cite-seq data
length(CITE_cells[CITE_cells %in% GEX_cells]) #how many from the gene expression?


pat4_SILP <- subset(pat4_SILP, cells = GEX_cells[GEX_cells %in% CITE_cells])
pat4_SILP[["CITE"]] <- CreateAssayObject(counts = pat04_SILP_CITE[,GEX_cells[GEX_cells %in% CITE_cells]])

#from QC cut offs
neg_object <- subset(pat4_SILP, subset = nFeature_RNA < 700)
pat4_SILP <- subset(pat4_SILP, subset = nFeature_RNA > 700 & nFeature_RNA < 5900 & percent.mt < 10)

# non sparse CITEseq data actually store better in a regular materix so the as.matrix() call is not memory intensive.
neg_adt_matrix = GetAssayData(neg_object, assay = "CITE", slot = 'counts') %>% as.matrix()
positive_adt_matrix = GetAssayData(pat4_SILP, assay = "CITE", slot = 'counts') %>% as.matrix()

isotypes <- c("Iso-IgG1","Iso-IgG2a","Iso-IgG2b")

normalized_matrix = DSBNormalizeProtein(cell_protein_matrix = positive_adt_matrix,
                                        empty_drop_matrix = neg_adt_matrix,
                                        use.isotype.control = TRUE,
                                        isotype.control.name.vec = isotypes)

colnames(normalized_matrix) <- str_replace(colnames(normalized_matrix),"-1","")
colnames(normalized_matrix) <- paste("P4_SILP_",colnames(normalized_matrix), sep="")
#save the normalized matrix to be able to upload to other objects
saveRDS(normalized_matrix,file= "Pat4SILPAllCells_CITEseq_dsbnorm.rds")

#### Adding concatenated DSB norm. CITE matrix for all samples with CITE-seq to Seurat object ####
MNP <- readRDS("/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R1_AllCells_MNP_LP.rds")
DSBmat_pat4SILP <- readRDS("Pat4SILPAllCells_CITEseq_dsbnorm.rds")
DSBmat_pat4cLP <- readRDS("Pat4SILPAllCells_CITEseq_dsbnorm.rds")
DSBmat_pat5cLP <- readRDS("Pat4SILPAllCells_CITEseq_dsbnorm.rds")
DSBmat_pat6cLP <- readRDS("Pat4SILPAllCells_CITEseq_dsbnorm.rds")

## the sample without isotypes and CD209 and CD11c, pat5 cLP
notCITE_mat_pat5 <- rep(NA, 5*length(colnames(DSBmat_pat5cLP)))
dim(notCITE_mat_pat5) <- c(5,length(colnames(DSBmat_pat5cLP)))
colnames(notCITE_mat_pat5) <- colnames(DSBmat_pat5cLP)
rownames(notCITE_mat_pat5) <- rownames(DSBmat_pat4cLP)[!rownames(DSBmat_pat4cLP) %in% rownames(DSBmat_pat5cLP)]
DSBmat_pat5cLP <- rbind(DSBmat_pat5cLP, notCITE_mat_pat5)

# co catenate the cite matrices
DSBmat <- DSBmat_pat4SILP %>% 
  cbind(DSBmat_pat4cLP) %>% 
  cbind(DSBmat_pat5cLP) %>%
  cbind(DSBmat_pat6cLP)
saveRDS(DSBmat, file="Citeseq_DSBNormMat_AllCells.rds")

Inc_cite <- rownames(MNP_LP@meta.data)
Not_NA <- citenames[citenames %in% rownames(MNP_LP@meta.data)]
Tobe_NA <- Inc_cite[!Inc_cite %in% Not_NA]

## Samples without cite-seq data
notCITE_mat <- rep(NA, length(rownames(DSBmat_pat4cLP))*length(Tobe_NA))
dim(notCITE_mat) <- c(length(rownames(normalized_matrix)),length(Tobe_NA))
colnames(notCITE_mat) <- Tobe_NA
rownames(notCITE_mat) <- rownames(DSBmat_pat4cLP)

## Bind both citeseq and NA data together and choose only cells existing in MNP_LP object
DSBmat_MNP <- cbind(DSBmat,notCITE_mat)
dim(DSBmat_MNP) #42506 cells
saveRDS(DSBmat_MNP, file="Citeseq_DSBNormMat_AllCellsInclNA.rds")
# subset to only include cells from MNP object
DSBmat_MNP <- DSBmat_MNP[,Inc_cite]
#check that there are the same number of cells in the objects
dim(DSBmat_MNP)
MNP_LP

## create CITE assay and add data
MNP_LP[["CITE"]] <- CreateAssayObject(data = DSBmat_MNP)
MNP_LP <- ScaleData(MNP_LP, assay = "CITE")

# add metadata to distinguish cells with CITE data and another layer for cells with the additional ABs added (CD209, CD11c, + isotype)
MNP_LP@meta.data$cite <- "NO"
MNP_LP@meta.data[Not_NA,]$cite <- "YES"
MNP_LP@meta.data$cite2 <- "NO"
MNP_LP@meta.data[Not_NA[!Not_NA %in% colnames(DSBmat_pat5cLP)],]$cite2 <- "YES"


