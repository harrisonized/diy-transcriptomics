# This follows DIY_scRNAseq.R. This script walks through the quality assessment 
# (QA) and analysis of single cell RNA-seq data

# In the 1st 1/2 of the script, we'll practice some basics
# using a small (~1000 cell) dataset from human peripheral
# blood mononuclear cells (PBMCs). This dataset comes from
# the public datasets on the 10X Genomics website:
# https://www.10xgenomics.com/resources/datasets

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
source(file.path(wd, 'R', 'utils.R'))
source(file.path(wd, 'R', 'functions.R'))

library('Seurat')
library('Matrix')
library('DropletUtils')

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
log <- log_open(paste0("covid19_scrnaseq_seurat"
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))



# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading data...'))

raw_mtx <- Matrix::readMM(
    file.path(wd, 'data', 'covid19_scrnaseq', 'pbmc_1k_v3_scRNAseq_processed',
              'counts_unfiltered', 'cellranger', 'matrix.mtx')
)

# add ensemble gene_ids to the data matrix as rownames
genes <- read.csv(
    file.path(wd, 'data', 'covid19_scrnaseq', 'pbmc_1k_v3_scRNAseq_processed',
              'counts_unfiltered', 'cellranger', 'genes.tsv'),
    sep = '\t', header = FALSE
)
rownames(raw_mtx) <- genes[,1]

# add cell barcodes as column names
barcodes <- read.csv(
    file.path(wd, 'data', 'covid19_scrnaseq', 'pbmc_1k_v3_scRNAseq_processed',
              'counts_unfiltered', 'cellranger', 'barcodes.tsv'),
    sep = '\t', header = FALSE
)
colnames(raw_mtx) <- barcodes[,1]


# ----------------------------------------------------------------------
# Filter empty drops

log_print(paste(Sys.time(), 'Filtering empty drops...'))

# use DropletUtils package to get probability that each barcode is a cell
# this takes a minute
out <- DropletUtils::emptyDrops(raw_mtx)
# set threshold probability for calling a cell
keep = ((out$FDR <= 0.05) & (is.na(out$FDR)==FALSE))
filt_mtx <- raw_mtx[,keep] 

# write out filtered results
if (save==TRUE) {
    write10xCounts(
        file.path(wd, 'data', 'covid19_scrnaseq', 'pbmc_1k_v3_scRNAseq_processed',
                  'counts_filtered', 'cellranger'),
        filt_mtx,
        gene.symbol = genes[,2],
        overwrite=T
    )
}


# ----------------------------------------------------------------------
# Generate QA report

# this report will contain some useful metrics as well as the traditional
# log-transformed UMI rank plot (a.k.a. 'waterfall' plot)
# plot was first described in the Drop-seq paper: - Macosko et al. 2015, DOI:10.1016/j.cell.2015.05.002
# this plot has two important points that we will try to identify:
# 1. 'knee point' - is the point where the signed curvature is minimized. 
# This corresponds to a transition between a distinct subset of barcodes
# with large totals and the majority of barcodes with smaller totals
# 2. 'inflection point' - is the point on the curve where the first derivative is minimized. 
# This corresponds to the point past which cells cannot reliably be distinguished from background

log_print(paste(Sys.time(), 'Generating QA report...'))

# load run info from JSON files produced by Kb
kb_stats <- c(
    fromJSON(file = file.path(
        wd, 'data', 'covid19_scrnaseq', 'pbmc_1k_v3_scRNAseq_processed',
        'inspect.json')
    ), 
    fromJSON(file = file.path(
        wd, 'data', 'covid19_scrnaseq', 'pbmc_1k_v3_scRNAseq_processed',
        'run_info.json')
    )
) 

# determine chemistry version
tech <- grep('10X(.*)', strsplit(kb_stats$call, '\\s')[[1]], value=T) 

# make a nice/simple table that summarizes that stats
seq_stats <- data.frame(
    # get sequencing/alignment stats 
    stat = c('Sequencing technology', 'Number of reads processed',
             '% reads pseudoaligned', '% reads on whitelist'),
    value = prettyNum(
        c(tech, kb_stats$n_processed, kb_stats$p_pseudoaligned, 
          round(kb_stats$percentageReadsOnWhitelist,2)),
    big.mark = ',')
)

# calculate cell stats and save to df
pct_counts_in_cells <- round((sum(filt_mtx)/sum(raw_mtx))*100, 2) 
med_counts_cell <- median(colSums(filt_mtx))
med_genes_cell <- median(apply(filt_mtx, 2, function(x) sum(x >= 1)))
tot_genes_detected <- sum(rowSums(filt_mtx)>=1)
cell_stats <- data.frame(
    stat = c('Estimated number of cells', '% counts in cells', 
             'Median counts per cell', 'Median genes per cell',
             'Total genes detected'), 
    value = prettyNum(
        c(ncol(filt_mtx), pct_counts_in_cells, med_counts_cell,
          med_genes_cell, tot_genes_detected),
    big.mark = ',')
)

# get rank stats
stats <- barcodeRanks(raw_mtx)

# load raw cells
raw_cells <- read.csv(
    file.path(wd, 'data', 'covid19_scrnaseq', 'pbmc_1k_v3_scRNAseq_processed',
              'counts_unfiltered', 'cellranger', 'barcodes.tsv'),
    sep = '\t', header = FALSE
)
filt_cells <- read.csv(
    file.path(wd, 'data', 'covid19_scrnaseq', 'pbmc_1k_v3_scRNAseq_processed',
              'counts_filtered', 'cellranger', 'barcodes.tsv'),
    sep = '\t', header = FALSE
)

# create barcode rank plot png
outfile = file.path(wd, 'figures', 'covid19_scrnaseq', 'barcode_rank.png')
if (save==TRUE) {
    if (!dir.exists(dirname(outfile))) {
      dir.create(dirname(outfile))
    }
    bc_rank_plot(
        stats = stats, raw_cells = raw_cells, filt_cells = filt_cells,
        save = outfile
    )
}


# output a HTML summary of the run
# currently broken because references in the function are off
# print_HTML(
#     seq_stats = seq_stats, cell_stats = cell_stats,
#     dir = file.path(wd, 'data', 'covid19_scrnaseq', 'pbmc_1k_v3_scRNAseq_processed',
#               'counts_unfiltered'),
#     sample_id = NULL
# )


# ----------------------------------------------------------------------
# Seurat

log_print(paste(Sys.time(), 'Creating Seurat object...'))

expression_matrix <- Read10X(
    file.path(wd, 'data', 'covid19_scrnaseq', 'pbmc_1k_v3_scRNAseq_processed',
              'counts_filtered', 'cellranger'),
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

# actually creating the Seurat Object
pbmc.1k.seurat <- CreateSeuratObject(counts = expression_matrix, min.cells = 3) %>% 
    NormalizeData(verbose = FALSE) %>% 
    FindVariableFeatures(verbose = FALSE)

# Let's calculate percent of mitochondrial reads
# NOTE: change 'MT' to 'mt' for mouse
pbmc.1k.seurat[["percent.mt"]] <- PercentageFeatureSet(object = pbmc.1k.seurat, pattern = "^MT-") 

VlnPlot(pbmc.1k.seurat, c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0.1)
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'covid19_scrnaseq', 'violin-pct_mitochondrial_reads.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# Filter your data
# NOTE: you need to be careful when setting cut-offs that you're not losing unique cell populations
pbmc.1k.seurat <- subset(
    pbmc.1k.seurat,
    subset = nCount_RNA < 20000 & 
             nCount_RNA > 1000 & 
             nFeature_RNA > 1000 & 
             percent.mt < 20
)

# another QA plot
# Potential things to look for in the type of QA plot produced above:
# 1. Data points in the bottom LEFT hand quadrant = low genes and UMIs per cell.
#    May represent poor quality cells.
# 2. Data points in the bottom RIGHT hand quadrant = low genes but high UMIs per cell.
#    These could be dying cells, but also could represent a population of a low complexity
#    celltype (i.e red blood cells).
ggplot(pbmc.1k.seurat@meta.data, aes(nCount_RNA, nFeature_RNA)) +
    geom_point(alpha = 0.7, size = 0.5) +
    labs(x = "Total UMI counts per cell", y = "Number of genes detected")
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'covid19_scrnaseq', 'scatter-num_genes_vs_counts.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Plot UMAP

log_print(paste(Sys.time(), 'UMAP Clustering...'))

# it is standard practice to apply a linear transformation ('scaling') before PCA.
# For single cell data this includes:
# 1. Shifting the expression of each gene, so that the mean expression across cells is 0
# 2. Scaling the expression of each gene, so that the variance across cells is 1
# This gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
pbmc.1k.seurat <- ScaleData(pbmc.1k.seurat, verbose = FALSE)
pbmc.1k.seurat <- RunPCA(pbmc.1k.seurat, npcs = 40, verbose = FALSE)
pbmc.1k.seurat <- RunUMAP(pbmc.1k.seurat, reduction = "pca", dims = 1:40)
pbmc.1k.seurat <- FindNeighbors(pbmc.1k.seurat, reduction = "pca", dims = 1:40)
pbmc.1k.seurat <- FindClusters(pbmc.1k.seurat, resolution = 0.5)

DimPlot(pbmc.1k.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'covid19_scrnaseq', 'umap.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Find cluster-specific genes

# generally speaking there are three main ways you can find cluster-specific marker genes with Seurat
# 1. 'FindMarkers' to compare a select cluster to all other cells not in that cluster
# 2. 'FindAllMarkers' to compare EACH cluster to all other cells
# 3. 'FindConservedMarkers' to identify genes conserved (shared) between two defined clusters

# We'll start with FindMarkers, since it allows you to choose exactly which cluster you'd like to focus on.
cluster1.markers <- FindMarkers(pbmc.1k.seurat, ident.1 = 1, min.pct = 0.25)
cluster1.markers$pct.diff <- cluster1.markers$pct.1 - cluster1.markers$pct.2
cluster1.markers.df <- as_tibble(cluster1.markers, rownames = "geneID")
# Export DEGs for each cluster (ranked by avg_logFC > 0.5)
myTopHits_cluster1 <- cluster1.markers.df %>% arrange(desc(avg_log2FC))
myTopHits_cluster1 <- dplyr::slice(myTopHits_cluster1, 1:20)

# interactive table
# datatable(myTopHits_cluster1, 
#           extensions = c('KeyTable', "FixedHeader"), 
#           caption = 'Table 1: Cluster 1 genes',
#           options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10,
#                          lengthMenu = c("10", "25", "50", "100"))) %>%
#     formatRound(columns=c(2:11), digits=2)

if (save==TRUE) {
    write_tsv(
        myTopHits_cluster1,
        file.path(wd, 'data', 'covid19_scrnaseq', 'pbmc_1k_v3_scRNAseq_processed',
                 'cluster1.tsv'),
    )
}

# plot genes of interest on UMAP
FeaturePlot(pbmc.1k.seurat, 
            reduction = "umap", 
            features = c("IGHM"),
            pt.size = 0.4, 
            order = TRUE,
            #split.by = "orig.ident",
            min.cutoff = 'q10',
            label = FALSE)
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'covid19_scrnaseq', 'feature_plot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# now let's try with FindAllMarkers
pbmc.1k.markers <- FindAllMarkers(pbmc.1k.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# let's take the top 10 marker genes for each cluster and plot as a heatmap
top10 <- pbmc.1k.markers %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc.1k.seurat, features = top10$gene)
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'covid19_scrnaseq', 'heatmap.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
