### Rscript for running imputation and tspace trajectory ####
## based on Line Wulff, 20.08.04
## Camilla Lemvigh, 13.12.2024

## libraries
library(Seurat)
library(magicBatch)
library(tSpace)
library(tidyverse)

dir <- "MNP_scRNAseq_paper"
control_df <- read.csv(paste0(dir, "/control_df.csv", sep = ""), header = T, sep = ";")
control_df <- control_df[,c(1,2)]
colnames(control_df) <- c("parameter","value")

#### read in Seurat data ####
obj_name <- control_df |> filter(parameter == "Seurat_obj") |> pull(value)
obj <- readRDS(paste0(dir, "/", obj_name, ".rds", sep = ""))
DefaultAssay(obj) <- "RNA"
obj <- FindVariableFeatures(obj, nfeatures = 1000)

# remove precursors
Idents(obj) <- "cDC_annotations"
obj <- subset(obj, idents = "proHLAlow", invert = T)

##### magic imputation on normalized data ####
impute <- control_df |> filter(parameter == "impute") |> pull(value)
if (impute == TRUE){
  #extract data from Seurat object
  #decide which genes are to be imputated
  if (control_df[control_df$parameter=="magic_features",]$value=="ALL"){
    var_feat <- rownames(obj[[obj@active.assay]])} else {
      var_feat <- VariableFeatures(obj)}
  
  DefaultAssay(obj) <- "RNA"
  magic_data <- FetchData(obj,vars = var_feat, slot="counts")
  magic_data <- as.matrix(magic_data)
  aff_magic_data <- obj@reductions$pca@cell.embeddings[,1:as.numeric(as.character(control_df[control_df$parameter=="pca_incl",]$value))]
  
  #magic imputation
  magic_data_out <- magicBatch(data=magic_data, aff_mat_input=aff_magic_data, t=6, #oversmoothed for tspace
                               python_command = "python3")
  #extract imputed matrix
  imputed <- magic_data_out[[1]]
  print("imputation done, 5x5 matrix:"); print(imputed[1:5,1:5])
  
  #add imputation data to Seurat obj
  #if (control_df[control_df$parameter=="save_imputation",]$value==TRUE){
  #magic_data <- magicBatch(data=magic_data, aff_mat_input=aff_magic_data, t=2, #save only t=2
  #                         python_command = "python3")
  #imputed <- magic_data[[1]]
  #obj[["imputated"]] <- CreateAssayObject(data = t(imputed))
  #saveRDS(obj, file = paste(dir,"/output/",obj_name,"_imputed.rds",sep=""))}
} #if impute == T

#### tSPACE ####
tspace <- control_df |> filter(parameter == "run_tspace") |> pull(value)
if (tspace == TRUE){
  #tspace_df <- FALSE
  #which input to use, imputed variable features or adjusted pca
  input <- control_df |> filter(parameter == "tspace_features") |> pull(value)
  
  if (input == "variable"){
    #obj <- FindVariableFeatures(obj, nfeatures = 3000)
    var_genes <- VariableFeatures(obj)
    var_genes <- var_genes[var_genes %in% colnames(imputed)]
    tspace_df <- imputed[, var_genes]
    print("The number of genes used for tspace is:");print(dim(tspace_df)[2])
    dimred <- "both"
  } else if (input == "pca"){
    pca_dims <- control_df |> filter(parameter == "pca_incl") |> pull(value) |> as.numeric()
    tspace_df <- obj@reductions$pca@cell.embeddings[, 1:pca_dims]
    dimred <- "both"
  } else if (input == "harmony"){
    pca_dims <- control_df |> filter(parameter == "pca_incl") |> pull(value) |> as.numeric()
    tspace_df <- obj@reductions$harmony@cell.embeddings[, 1:pca_dims]
    dimred <- "both"
  } else {
    print("You have not assigned a matrix for the tSpace calculation")
  }
  #if (tspace_df!=FALSE){
  class(tspace_df)
  tspace_out <- tSpace(tspace_df, K = 20, graph = 5,
                       trajectories = 100, wp = 20, ground_truth = F,
                       weights = "exponential", dr = dimred, seed = 8, core_no = 16)
  #save files
  print(paste0(dir, "/tspace/", obj_name, "_integrated_sub_tspacefile.rds",sep=""))
  #tspace_out
  saveRDS(tspace_out, file = paste0(dir, "/tspace/", obj_name, "_integrated_sub_tspacefile.rds", sep = ""))
  #tspace_out2 <- tspace_out$ts_file
  #saveRDS(tspace_out2,file = paste(dir,"/output/",obj_name,"_tspace_tsfile.rds",sep=""))
  #} #tspace_df cannot be empty
} #end if run_tspace is true

print("All done!")



