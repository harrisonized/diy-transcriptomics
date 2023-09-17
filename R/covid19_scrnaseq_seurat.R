## Adapted from: DIY_scRNAseq.R
##
## In the 1st 1/2 of the script, we'll practice some basics
## using a small (~1000 cell) dataset from human peripheral
## blood mononuclear cells (PBMCs). This dataset comes from
## the public datasets on the 10X Genomics website:
## https://www.10xgenomics.com/resources/datasets

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
source(file.path(wd, 'R', 'utils.R'))
source(file.path(wd, 'R', 'functions.R'))  # 
suppressMessages(library('Seurat'))
library('Matrix')
suppressMessages(library('DropletUtils'))
library('tibble')
suppressMessages(library('dplyr'))
library('rjson')
library('R2HTML')
library('readr')
suppressMessages(library('textshape'))
library('scales')
library('ggplot2')
suppressMessages(library('DT'))
library('optparse')
library('logr')


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"),
                default='data/covid19_scrnaseq/pbmc_1k_v3_scRNAseq_processed',
                metavar='data/covid19_scrnaseq/pbmc_1k_v3_scRNAseq_processed',
                type="character",
                help="cellranger directory containing barcodes.tsv, genes.tsv, and matrix.mtx"),
 
    make_option(c("-o", "--output-dir"), default="figures/covid19_scrnaseq",
                metavar="figures/covid19_scrnaseq", type="character",
                help="set the output directory for the figures"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt['troubleshooting'][[1]]

if (!troubleshooting) {
    if (!dir.exists(file.path(wd, opt['output-dir'][[1]]))) {
        dir.create(file.path(wd, opt['output-dir'][[1]]))
    }
}

# Start Log
start_time = Sys.time()
log <- log_open(paste0("covid19_scrnaseq_seurat",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading data...'))

filt_mtx <- Read10X(
    file.path(wd, opt['input-dir'][[1]], 'counts_filtered', 'cellranger'),
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

# ----------------------------------------------------------------------
# Seurat

log_print(paste(Sys.time(), 'Creating Seurat object...'))

# actually create the Seurat Object
seurat_obj <- CreateSeuratObject(counts = filt_mtx, min.cells = 3) %>% 
    NormalizeData(verbose = FALSE) %>% 
    FindVariableFeatures(verbose = FALSE)

# Let's calculate percent of mitochondrial reads
# NOTE: change 'MT' to 'mt' for mouse
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-") 

VlnPlot(seurat_obj, c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0.1)
if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'violin-pct_mitochondrial_reads.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# Filter your data
# NOTE: you need to be careful when setting cut-offs that you're not losing unique cell populations
seurat_obj <- subset(
    seurat_obj,
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
ggplot(seurat_obj@meta.data, aes(nCount_RNA, nFeature_RNA)) +
    geom_point(alpha = 0.7, size = 0.5) +
    labs(x = "Total UMI counts per cell", y = "Number of genes detected")
if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'scatter-num_genes_vs_counts.png'),
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
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 40, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:40)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:40)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

DimPlot(seurat_obj, reduction = "umap", split.by = "orig.ident", label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'umap.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Find cluster-specific genes

# generally speaking there are three main ways you can find cluster-specific marker genes with Seurat
# 1. 'FindMarkers' to compare a select cluster to all other cells not in that cluster
# 2. 'FindAllMarkers' to compare EACH cluster to all other cells
# 3. 'FindConservedMarkers' to identify genes conserved (shared) between two defined clusters

# We'll start with FindMarkers, since it allows you to choose exactly which cluster you'd like to focus on.
cluster1.markers <- FindMarkers(seurat_obj, ident.1 = 1, min.pct = 0.25)
cluster1.markers[['pct.diff']] <- cluster1.markers[['pct.1']] - cluster1.markers[['pct.2']]
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

if (!troubleshooting) {
    write_tsv(myTopHits_cluster1, 
              file.path(wd, opt['input-dir'][[1]], 'cluster1.tsv'))
}

# plot genes of interest on UMAP
FeaturePlot(seurat_obj, 
            reduction = "umap", 
            features = c("IGHM"),
            pt.size = 0.4, 
            order = TRUE,
            #split.by = "orig.ident",
            min.cutoff = 'q10',
            label = FALSE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'feature_plot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# now let's try with FindAllMarkers
pbmc.1k.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# let's take the top 10 marker genes for each cluster and plot as a heatmap
top10 <- pbmc.1k.markers %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat_obj, features = top10[['gene']])
if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'heatmap.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
