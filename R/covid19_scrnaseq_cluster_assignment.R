# This follows DIY_scRNAseq.R. This script walks through the quality assessment 
# (QA) and analysis of single cell RNA-seq data

# required for cellassign to access tensorflow through python
reticulate::use_condaenv('r-reticulate')

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
library('DropletUtils')
library('scran')
library('cellassign')  # takes a few seconds
library('celldex') #a large collection of reference expression datasets with curated cell type labels for use with SingleR package
library('SingleR') #automated cell type annotation ('label transfer') using reference data

library('pheatmap')
# library('scater')  # plotUMAP, but could not install this

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


# ----------------------------------------------------------------------
# Use cellassign to assign clusters

# compute size factors
pbmc.1k.sce <- scran::computeSumFactors(pbmc.1k.sce)

fit <- cellassign(
  exprs_obj = sce_marker,
  marker_gene_info = pbmc_marker_matrix,
  s = sizeFactors(pbmc.1k.sce),
  shrinkage = TRUE,
  max_iter_adam = 50,
  min_delta = 2,
  verbose = TRUE)

# > fit
# A cellassign fit for 1223 cells, 22 genes, 5 cell types with 0 covariates
#             To access cell types, call celltypes(x)
#             To access cell type probabilities, call cellprobs(x)

# incorporate the cellAssign result into your singleCellExperiment
pbmc.1k.sce$cell_type <- fit$cell_type

# plotUMAP is the Scater equivalent of Seurat's DimPlot
# plotUMAP(pbmc.1k.sce, colour_by = "cell_type")  # do we really need this?


# ----------------------------------------------------------------------
# Use SingleR to assign clusters

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

plotScoreHeatmap(predictions)  # SingleR

# now add back to singleCellExperiment object (or Seurat objects)
pbmc.1k.sce[["SingleR.labels"]] <- predictions$labels
# plotUMAP(pbmc.1k.sce, colour_by = "SingleR.labels")


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
