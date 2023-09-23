## Adapted from: DIY_scRNAseq.R
## Creates a Seurat object from the empty-drop-filtered cellranger output.
## Plots UMAPs and heatmaps of all DEGs and cluster-specific DEGs

wd = dirname(this.path::here())  # wd = '~/github/R/diy-transcriptomics'
suppressMessages(library('Seurat'))
library('Matrix')
suppressMessages(library('DropletUtils'))
library('tibble')
suppressMessages(library('dplyr'))
library('rjson')
library('R2HTML')
library('readr')
suppressMessages(library('textshape'))
suppressMessages(library('scales'))
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

    make_option(c("-g", "--gene-of-interest"), default="IGHM",
                metavar="IGHM", type="character",
                help="choose a gene"),

    make_option(c("-c", "--cluster-id"), default=1, metavar="1",
                type="integer", help="set the cluster number for the cluster-specific analysis"),

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

seurat_obj <- CreateSeuratObject(counts = filt_mtx, min.cells = 3) %>% 
    NormalizeData(verbose = FALSE) %>% 
    FindVariableFeatures(verbose = FALSE)

# calculate percent mitochondrial reads
seurat_obj[["pct_mt_reads"]] <- PercentageFeatureSet(
    object = seurat_obj,
    pattern = "^MT-"  # NOTE: change 'MT' to 'mt' for mouse
)

# QC plot
VlnPlot(seurat_obj, c("nCount_RNA", "nFeature_RNA", "pct_mt_reads"), pt.size = 0.1)
if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'violin-raw_qc.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Filter data

# NOTE: you need to be careful when setting cut-offs that you're not losing unique cell populations
# Note: seurat_obj is an instance of a class, where attributes are accessed via the "@" operator
# See: http://adv-r.had.co.nz/S4.html
seurat_obj <- subset(
    seurat_obj,
    subset = ((nCount_RNA < 20000) & 
              (nCount_RNA > 1000) & 
              (nFeature_RNA > 1000) & 
              (pct_mt_reads < 20))
)

# another QC plot
ggplot(seurat_obj@meta.data, aes(nCount_RNA, nFeature_RNA)) +
    geom_point(alpha = 0.7, size = 0.5) +
    labs(x = "Total UMI counts per cell", y = "Number of genes detected")
# Potential things to look for in the type of QA plot produced above:
# 1. Data points in the bottom LEFT hand quadrant = low genes and UMIs per cell.
#    May represent poor quality cells.
# 2. Data points in the bottom RIGHT hand quadrant = low genes but high UMIs per cell.
#    These could be dying cells, but also could represent a population of a low complexity
#    celltype (i.e red blood cells).

if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'scatter-num_genes_vs_counts.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Standard Clustering Workflow

log_print(paste(Sys.time(), 'Running standard clustering workflow...'))

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

# Plot UMAP with clusters highlighted
DimPlot(seurat_obj, reduction = "umap", split.by = "orig.ident", label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'umap-clusters.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# Plot UMAP for specific gene of interest
gene = opt['gene-of-interest'][[1]]  # gene = 'IGHM'
FeaturePlot(
    seurat_obj, reduction = "umap", features = c(gene),
    pt.size = 0.4,  min.cutoff = 'q10',
    # split.by = "orig.ident",
    order = TRUE, label = FALSE
)
if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], paste0('umap-', gene, '.png')),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Find cluster-specific DEGs

cluster_id = opt['cluster-id'][[1]]  # cluster_id = 1

# generally speaking there are three main ways you can find cluster-specific marker genes with Seurat
# 1. 'FindMarkers' to compare a select cluster to all other cells not in that cluster
# 2. 'FindAllMarkers' to compare EACH cluster to all other cells
# 3. 'FindConservedMarkers' to identify genes conserved (shared) between two defined clusters
# See: https://satijalab.org/seurat/reference/findmarkers

cluster_degs <- FindMarkers(seurat_obj, ident.1 = cluster_id, min.pct = 0.25)
cluster_degs[['pct.diff']] <- cluster_degs[['pct.1']] - cluster_degs[['pct.2']]
cluster_degs <- as_tibble(cluster_degs, rownames = "geneID")

# Export top 20 DEGs for cluster 1
cluster_degs_top20 <- cluster_degs %>% arrange(desc(avg_log2FC)) %>% slice(1:20)
if (!troubleshooting) {
    write_tsv(
        cluster_degs_top20,
        file.path(wd, opt['input-dir'][[1]],
                  paste0('degs-top_20-cluster_', opt['cluster-id'][[1]], '.tsv'))
    )
}

# # interactive table
# datatable(degs_top20,
#           extensions = c('KeyTable', "FixedHeader"), 
#           caption = paste('Table 1: Cluster', opt['cluster-id'][[1]] , 'genes'),
#           options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10,
#                          lengthMenu = c("10", "25", "50", "100"))
#     ) %>% formatRound(columns=c(2:11), digits=2)


# ----------------------------------------------------------------------
# Plot heatmap for all DEGs

all_degs <- FindAllMarkers(seurat_obj,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Export top 10 DEGs for each cluster and plot as a heatmap
all_degs_top10 <- all_degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat_obj, features = all_degs_top10[['gene']])
if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'heatmap.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
