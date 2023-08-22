# This follows DIY_scRNAseq.R. This script walks through the quality assessment 
# (QA) and analysis of single cell RNA-seq data

# In the 2nd 1/2 of the script, we'll import two separate Seurat
# objects generated from the spleen of naive and Toxoplasma gondii
# infected mice, giving us an opportunity to create and analyze an
# integrated dataset

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
source(file.path(wd, 'R', 'utils.R'))
source(file.path(wd, 'R', 'functions.R'))

library('Seurat')
library('Matrix')
library('DropletUtils')

reticulate::use_condaenv('r-reticulate')  # required for cellassign to access tensorflow through python
library('cellassign')  # takes a few seconds

library('scran')
library('SingleR') #automated cell type annotation ('label transfer') using reference data
library('celldex') #a large collection of reference expression datasets with curated cell type labels for use with SingleR package

# library('scater')  # could not install  # what do we need this for?

# library(tidyverse)  # bad practice
library('tibble')
library('dplyr')
library('rjson')
library('R2HTML')
library('readr')
library('textshape')
library('scales')

library('ggplot2')
library('DT')
library('pheatmap')

library('optparse')
library('logr')


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-s", "--save"), default=TRUE, action="store_false", metavar="TRUE",
                type="logical", help="disable if you're troubleshooting and don't want to overwrite your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

save = opt['save'][[1]]  # save=FALSE

# Start Log
start_time = Sys.time()
log <- log_open(paste("eda ", start_time, '.log', sep=''))
log_print(paste('Script started at:', start_time))



# ----------------------------------------------------------------------
# Read Data

# Note: to convert Seurat object to SingleCellExperiment object directly:
# pbmc.1k.sce <- as.SingleCellExperiment(pbmc.1k.seurat)

log_print(paste(Sys.time(), 'Reading data...'))
pbmc.1k.sce <- read10xCounts(
    file.path(wd, 'data', 'covid19_scrnaseq', 'pbmc_1k_v3_scRNAseq_processed',
          'counts_filtered', 'cellranger')
)

 # have row names be gene_name instead of gene_id
rownames(pbmc.1k.sce) <- rowData(pbmc.1k.sce)$Symbol

# inspect
# colnames(pbmc.1k.sce)  # null
# reducedDims(pbmc.1k.sce)  # null
# assays(pbmc.1k.sce)  # counts


# ----------------------------------------------------------------------
# Assign identity to cell clusters

# create a list of markers
# you can find cell specific markers here: http://biocc.hrbmu.edu.cn/CellMarker/
pbmc_marker_list <- list(
  `Monocytes` = c("CD14", "CD68"),
  `T cells` = c("CD2", "CD3D", "TRAC", "IL32", "CD3E", "PTPRC"),
  `NK cells` = c("GZMK", "KLRF1", "CCL3", "CMC1", "NKG7", "PTPRC"),
  `Plasma cells` = c("CD27", "IGHG1", "CD79A", "IGHG2", "PTPRC", "IGKC"),
  `Mature B cells` = c("MS4A1", "LTB", "CD52", "IGHD", "CD79A", "PTPRC", "IGKC")
)

# convert your marker gene list from above to a matrix
pbmc_marker_matrix <- marker_list_to_mat(pbmc_marker_list, include_other = FALSE)

# you can view this matrix as a heatmap
if (save==TRUE) {
    filename=file.path(wd, 'figures', 'covid19_scrnaseq', 'pheatmap_example_1.png')
} else {
    filename=NA
}
pheatmap(pbmc_marker_matrix,
         filename=filename  # save
)


# ----------------------------------------------------------------------
# Subset data

# subset
my.subset <- pbmc.1k.sce[,c(1,2,8)]

# make sure all your markers were actually observed in your single cell data.  Remove markers that were not detected
marker_in_sce <- match(rownames(pbmc_marker_matrix), rowData(pbmc.1k.sce)$Symbol)
# stopifnot(all(!is.na(marker_in_sce)))

# subset data to include only markers
sce_marker <- pbmc.1k.sce[marker_in_sce, ]
# stopifnot(all.equal(rownames(pbmc_marker_matrix), rowData(sce_marker)$Symbol))

# compute size factors
pbmc.1k.sce <- scran::computeSumFactors(pbmc.1k.sce)

# run cellAssign
fit <- cellassign(
  exprs_obj = sce_marker,
  marker_gene_info = pbmc_marker_matrix,
  s = sizeFactors(pbmc.1k.sce),
  shrinkage = TRUE,
  max_iter_adam = 50,
  min_delta = 2,
  verbose = TRUE)

quit()



# incorporate the cellAssign result into your singleCellExperiment
pbmc.1k.sce$cell_type <- fit$cell_type
# plotUMAP is the Scater equivalent of Seurat's DimPlot
plotUMAP(pbmc.1k.sce, colour_by = "cell_type")

# a different way of labeling clusters using public datasets
# now label using singleR and celldex (requires an internet connection to connect to ExperimentHub)
ENCODE.data <- BlueprintEncodeData(ensembl = FALSE) #259 RNA-seq samples of pure stroma and immune cells as generated and supplied by Blueprint and ENCODE
HPCA.data <- HumanPrimaryCellAtlasData(ensembl = FALSE) #713 microarray samples from the Human Primary Cell Atlas (HPCA) (Mabbott et al., 2013).
DICE.data <- DatabaseImmuneCellExpressionData(ensembl = FALSE) #1561 bulk RNA-seq samples of sorted immune cell populations
ImmGen.data <- ImmGenData(ensembl = FALSE) # 830 microarray samples of pure mouse immune cells, generated by the Immunologic Genome Project (ImmGen)
Monaco.data <- MonacoImmuneData(ensembl = FALSE) #114 bulk RNA-seq samples of sorted immune cell populations that can be found in GSE107011.
MouseRNAseq.data <- MouseRNAseqData(ensembl = FALSE) #358 bulk RNA-seq samples of sorted cell populations
Hemato.data <- NovershternHematopoieticData(ensembl = FALSE) #211 bulk human microarray samples of sorted hematopoietic cell populations that can be found in GSE24759


predictions <- SingleR(test=pbmc.1k.sce, assay.type.test=1, 
                       ref=Monaco.data, labels=Monaco.data$label.main)

plotScoreHeatmap(predictions)

#now add back to singleCellExperiment object (or Seurat objects)
pbmc.1k.sce[["SingleR.labels"]] <- predictions$labels
plotUMAP(pbmc.1k.sce, colour_by = "SingleR.labels")

# Integrate multiple scRNA-seq datasets ----
# To demonstrate integration, we'll leave behind the PBMC dataset we worked with above
# We'll read in two Seurat objects - one generated from the spleen of a untreated mouse (control), and the second from the spleen of mouse infected with Toxoplasma gondii
load("spleen.naive.seurat")
DimPlot(spleen.naive.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)

load("spleen.toxoInfected.seurat")
DimPlot(spleen.toxoInfected.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)

# since we are now going to work with multiple samples, we need a study design file with our sample metadata
targets <- read_tsv("studyDesign.txt")

# extract variables of interest
sampleID <- targets$sampleID
treatment <- targets$treatment

# annotate your seurat objects with as much or as little metadata as you want!
spleen.naive.seurat$treatment <- treatment[1]
spleen.toxoInfected.seurat$treatment <- treatment[2]

# take a look at where this metadata lives in the seurat object
spleen.toxoInfected.seurat@meta.data$treatment

# select features that are repeatedly variable across datasets for integration
spleen_features <- SelectIntegrationFeatures(object.list = c(spleen.naive.seurat, spleen.toxoInfected.seurat))
spleen_anchors <- FindIntegrationAnchors(object.list = c(spleen.naive.seurat, spleen.toxoInfected.seurat), anchor.features = spleen_features)
spleen_integrated <- IntegrateData(anchorset = spleen_anchors)
# NOTE: if you look at your seurat object, the default assay as changed from 'RNA' to 'integrated'
# this can be change anytime using the line below
# this would be the same way you would change between scRNA-seq and scATAC-seq
# DefaultAssay(spleen_integrated) <- "RNA"

# Run the standard workflow for visualization and clustering
spleen_integrated <- ScaleData(spleen_integrated, verbose = FALSE)
spleen_integrated <- RunPCA(spleen_integrated, npcs = 30, verbose = FALSE)
spleen_integrated <- RunUMAP(spleen_integrated, reduction = "pca", dims = 1:30)
spleen_integrated <- FindNeighbors(spleen_integrated, reduction = "pca", dims = 1:30)
spleen_integrated <- FindClusters(spleen_integrated, resolution = 0.5)
DimPlot(spleen_integrated, reduction = "umap", label = TRUE)

# let's see what proportion of our total cells reside in each cluster
prop.table(table(Idents(spleen_integrated)))

# remember, we have metadata in this integrated seurat object, so you can use this to split your UMAP
DimPlot(spleen_integrated, reduction = "umap", 
        split.by = "treatment", # this facets the plot 
        group.by = "seurat_clusters", # labels the cells with values from your group.by variable
        label = TRUE)

# plot genes of interest on UMAP
FeaturePlot(spleen_integrated, 
            reduction = "umap", 
            features = 'Sdc1',
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE)

# we can plot more than one gene here
my_fav_genes <- c("Cd4", "Cd8a")
FeaturePlot(spleen_integrated, 
            reduction = "umap", 
            features = my_fav_genes,
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE)

# Leveraging cluster identity in your analysis ----
# now let's rerun our cluster identification using SingleR
spleen_integrated.sce <- as.SingleCellExperiment(spleen_integrated)
predictions <- SingleR(test=spleen_integrated.sce, assay.type.test=1, 
                       ref=MouseRNAseq.data, labels=MouseRNAseq.data$label.main)

#now add back to singleCellExperiment object (or Seurat objects)
spleen_integrated.sce[["SingleR.labels"]] <- predictions$labels
plotUMAP(spleen_integrated.sce, colour_by = "SingleR.labels")

spleen_integrated2 <- as.Seurat(spleen_integrated.sce, counts = NULL)
DimPlot(spleen_integrated2, reduction = "UMAP", 
        split.by = "treatment", # this facets the plot 
        group.by = "SingleR.labels", # labels the cells with values from your group.by variable
        label = TRUE)

# If we repeat the steps above for different cell type markers, we get a sense for the following cluster IDs
new.cluster.ids <- c("B cells", "RBCs", "CD8+ T cells", "B cells", "RBCs", "CD4+ T cells", "CD4+ T cells", "Monocytes/Macrophages", "Granulocytes", "Monocytes/Macrophages", "B cells", "Plasma cells", "Monocytes/Macrophages", "Monocytes/Macrophages", "Granulocytes", "CD8+ T cells", "CD8+ T cells", "17", "18", "19", "20") 
names(new.cluster.ids) <- levels(spleen_integrated)
spleen_integrated <- RenameIdents(spleen_integrated, new.cluster.ids)
DimPlot(spleen_integrated, reduction = "umap", 
        split.by = "treatment", # this facets the plot 
        label = TRUE)

# take a look at what you've done
Idents(spleen_integrated)

# subset seurat object to focus on single cluster ----
# let's get just the CD4 T cells
spleen_integrated.CD4.Tcells <- subset(spleen_integrated, idents = "CD4+ T cells")
# you could re-cluster if you want...depends on what you're trying to accomplish
#spleen_integrated.CD4.Tcells <- FindNeighbors(spleen_integrated.CD4.Tcells, dims = 1:10, k.param = 5)
#spleen_integrated.CD4.Tcells <- FindClusters(spleen_integrated.CD4.Tcells)
DimPlot(spleen_integrated.CD4.Tcells, reduction = "umap", label = TRUE)
Idents(spleen_integrated.CD4.Tcells)

# now we need to switch out 'Idents' to be treatment, rather than cluster
Idents(spleen_integrated.CD4.Tcells) <- spleen_integrated.CD4.Tcells$treatment
inf.vs.naive.markers <- FindMarkers(object = spleen_integrated.CD4.Tcells, 
                                    ident.1 = "infected", 
                                    ident.2 = "naive", 
                                    min.pct = 0)

inf.vs.naive.markers$pct.diff <- inf.vs.naive.markers$pct.1 - inf.vs.naive.markers$pct.2
inf.vs.naive.markers.df <- as_tibble(inf.vs.naive.markers, rownames = "geneID")
# Export DEGs for each cluster (ranked by avg_logFC > 0.5)
myTopHits <- inf.vs.naive.markers.df %>% arrange(desc(avg_log2FC))

FeaturePlot(spleen_integrated.CD4.Tcells, 
            reduction = "umap", 
            features = "Ccl5",
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE)



end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
