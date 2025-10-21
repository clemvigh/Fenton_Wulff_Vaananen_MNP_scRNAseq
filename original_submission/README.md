# FentonWulff_LP_MNP
This repository contains all the most essential scripts, which were used in analysing the single cell data for the preprint **Fenton, Wulff et al. BioRxiv 2021**.

QC of single cells was first assessed per sample based on number of read counts, feature/gene counts and mitochondrial gene percentage. Barcodes determined to come from debris or doublets by these parameters were removed before the samples were integrated with Seurat's anchor integration and scaled while removing effects of cell cycle, MT gene load, read and gene depth. Based on 15 PC's UMAP and clustering were calculated and contaminating cell types were removed. The remaining cells were rescaled, pca rerun, reclustered and a new UMAP was run (still on 15 PCs).
The supercluster of CD14+ and CD1C+ cells was again subsetted, rescaled, reclustered etc. at resolution 2.8.

**CiteseqPerSample_AddCITEseqToMNPobj.R** Normalisation by denoising and scaling to background of CITE-seq data per sample. Adding CITE-seq data to the Seurat object.

**Running_tspace.R** - tSpace setup, requires a control_df.csv file, run on computerome (Danish HPC system) with anaconda3/4.0.0, R/4.0.0 Seurat/3.2.0. Examples of control_df.csv files for dendritic cell (Figure 2-4) and mono mac (Figure 5-6) analysis here named respectively DC_control_df.csv and monomac_control_df.csv.

**ProgenyDorotheaAnalysis.R** - example of how analysis with Dorothea and Progeny were run.

**GOplot_cDCexample_SupFig2.R** - Example of how the radarplot for GO analysis were created based on table outputs from EnrichR. All significant GO terms were manually curated to include immunologically relevant terms in the final figures.

**Fig4_putprecursors.R** - Analysis for Figure 4 and Supplementary Figure 4. csv files used in this script velocyto_cells_v3.csv, DC3prec_cl50.csv and DC1prec_cl5035.csv.
