## Adapted from: DIY_scRNAseq.R
## Assigns clusters and uses data from celldex to label them

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
reticulate::use_condaenv('r-reticulate')  # required for cellassign to access tensorflow through python
library('rjson')
library('cellassign')  # takes a few seconds
suppressMessages(library('DropletUtils'))
library('SingleR') # automated cell type annotation ('label transfer') using reference data
suppressMessages(library('celldex')) # a large collection of reference expression datasets with curated cell type labels for use with SingleR package
library('scran')
library('pheatmap')
# library('scater')  # plotUMAP, could not install this
library('optparse')
library('logr')


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"),
                default='data/covid19_scrnaseq/pbmc_1k_v3_scRNAseq_processed/counts_filtered/cellranger',
                metavar='data/covid19_scrnaseq/pbmc_1k_v3_scRNAseq_processed/counts_filtered/cellranger',
                type="character",
                help="cellranger directory containing barcodes.tsv, genes.tsv, and matrix.mtx"),

    make_option(c("-m", "--marker-list"),
                default='data/covid19_scrnaseq/pbmc_marker_list.json',
                metavar='data/covid19_scrnaseq/pbmc_marker_list.json',
                type="character",
                help="path/to/marker_list.json"),

    make_option(c("-c", "--celldex"),
                default='Monaco', metavar='Monaco', type="character",
                help="Choose from: ['ENCODE', 'HPCA', 'DICE', 'ImmGen', 'Monaco', 'MouseRNAseq', 'Hemato']"),

    make_option(c("-e", "--ensembl"),  default=FALSE,
                metavar="FALSE", action="store_true", type="logical",
                help="When querying celldex"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt['troubleshooting'][[1]]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("covid19_scrnaseq_cluster_assignment-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading data...'))

sce_counts <- read10xCounts(file.path(wd, opt['input-dir'][[1]]))
# Note: to convert Seurat object to SingleCellExperiment object directly:
# sce_counts <- as.SingleCellExperiment(pbmc.1k.seurat)
rownames(sce_counts) <- rowData(sce_counts)[['Symbol']]   # row names be gene_name instead of gene_id

marker_list <- fromJSON(file=file.path(wd, opt['marker-list'][[1]]))
marker_matrix <- marker_list_to_mat(marker_list, include_other = FALSE)  # convert to matrix

# view matrix as a heatmap
if (!troubleshooting) {
    filename=file.path(wd, 'figures', 'covid19_scrnaseq', 'pheatmap_example_1.png')
} else {
    filename=NA
}
pheatmap(
    marker_matrix,
    filename=filename  # save
)


# ----------------------------------------------------------------------
# Label clusters using cellassign

log_print(paste(Sys.time(), 'Labeling clusters using cellassign...'))

# only include markers detected in sce_counts
markers_in_sce <- match(rownames(marker_matrix), rowData(sce_counts)[['Symbol']])
sce_counts_subset <- sce_counts[markers_in_sce, ]  # filter sce_counts using markers
factors <- scran::computeSumFactors(sce_counts)  # compute size factors
fit <- cellassign(
    exprs_obj = sce_counts_subset,
    marker_gene_info = marker_matrix,
    s = sizeFactors(factors),
    shrinkage = TRUE,
    max_iter_adam = 50,
    min_delta = 2,
    verbose = TRUE
)
# A cellassign fit for 1223 cells, 22 genes, 5 cell types with 0 covariates
#             To access cell types, call celltypes(x)
#             To access cell type probabilities, call cellprobs(x)

# incorporate the cellAssign result into your singleCellExperiment
sce_counts[['cell_type']] <- fit[['cell_type']]

# Could not install scater
# plotUMAP(sce_counts, colour_by = "cell_type")


# ----------------------------------------------------------------------
# Label clusters using public datasets from celldex
# requires an internet connection to connect to ExperimentHub

# switch
switch=list(

    #259 RNA-seq samples of pure stroma and immune cells as generated and supplied by Blueprint and ENCODE
    'ENCODE'=BlueprintEncodeData,

    #713 microarray samples from the Human Primary Cell Atlas (HPCA) (Mabbott et al., 2013).
    'HPCA'=HumanPrimaryCellAtlasData,

    #1561 bulk RNA-seq samples of sorted immune cell populations
    'DICE'=DatabaseImmuneCellExpressionData,

    # 830 microarray samples of pure mouse immune cells, generated by the Immunologic Genome Project (ImmGen)
    'ImmGen'=ImmGenData,

    #114 bulk RNA-seq samples of sorted immune cell populations that can be found in GSE107011.
    'Monaco'=MonacoImmuneData,

    #358 bulk RNA-seq samples of sorted cell populations
    'MouseRNAseq'=MouseRNAseqData,

    #211 bulk human microarray samples of sorted hematopoietic cell populations that can be found in GSE24759
    'Hemato'=NovershternHematopoieticData
)
get_data=switch[opt['celldex'][[1]]][[1]]

ref_data <- get_data(ensembl = opt['ensembl'][[1]])  # ensembl=FALSE by default
predictions <- SingleR(
    test=sce_counts,
    assay.type.test=1,
    ref=ref_data,
    labels=ref_data[['label.main']]
)
plotScoreHeatmap(predictions)  # SingleR

# now add back to singleCellExperiment object (or Seurat objects)
sce_counts[["SingleR.labels"]] <- predictions[['labels']]
# plotUMAP(sce_counts, colour_by = "SingleR.labels")


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()


# ----------------------------------------------------------------------
# Code graveyard

# # create a list of markers, then export to JSON
# # you can find cell specific markers here: http://biocc.hrbmu.edu.cn/CellMarker/
# marker_list <- list(
#   "Monocytes" = c("CD14", "CD68"),
#   "T cells" = c("CD2", "CD3D", "TRAC", "IL32", "CD3E", "PTPRC"),
#   "NK cells" = c("GZMK", "KLRF1", "CCL3", "CMC1", "NKG7", "PTPRC"),
#   "Plasma cells" = c("CD27", "IGHG1", "CD79A", "IGHG2", "PTPRC", "IGKC"),
#   "Mature B cells" = c("MS4A1", "LTB", "CD52", "IGHD", "CD79A", "PTPRC", "IGKC")
# )

# if (!troubleshooting){
#     write(toJSON(marker_list), opt['marker-list'][[1]])
# }
