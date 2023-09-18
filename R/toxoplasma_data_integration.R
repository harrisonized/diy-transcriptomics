## Adapted from: DIY_scRNAseq.R
## Compare two Seurat objects from the spleen of naive mice versus
## Toxoplasma gondii infected mice.

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
suppressMessages(library('Seurat'))
library('Matrix')
suppressMessages(library('DropletUtils'))
reticulate::use_condaenv('r-reticulate')  # required for cellassign to access tensorflow through python
library('cellassign')  # takes a few seconds
library('scran')
library('SingleR') # automated cell type annotation ('label transfer') using reference data
library('celldex') # a large collection of reference expression datasets with curated cell type labels for use with SingleR package
# library('scater')  # plotUMAP, could not install this
# library(tidyverse)  # too broad
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
    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt['troubleshooting'][[1]]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("toxoplasma_data_integration-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Look at baseline data

log_print(paste(Sys.time(), 'Reading data...'))

# not a fan of this scheme
load(file.path(wd, 'data', 'toxoplasma', 'spleen.naive.seurat'))
DimPlot(spleen.naive.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'toxoplasma', 'umap_naive.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

load(file.path(wd, 'data', 'toxoplasma', 'spleen.toxoInfected.seurat'))
DimPlot(spleen.toxoInfected.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'toxoplasma', 'umap_infected.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Add annotations

log_print(paste(Sys.time(), 'Adding annotations...'))

# since we are now going to work with multiple samples, we need a study design file with our sample metadata
targets <- read_tsv(file.path(wd, 'data', 'toxoplasma', 'studyDesign.txt'))

# extract variables of interest
sampleID <- targets$sampleID
treatment <- targets$treatment

# annotate your seurat objects with as much or as little metadata as you want!
spleen.naive.seurat$treatment <- treatment[1]
spleen.toxoInfected.seurat$treatment <- treatment[2]


# ----------------------------------------------------------------------
# Integrate the data

# take a look at where this metadata lives in the seurat object
# spleen.toxoInfected.seurat@meta.data$treatment

# select features that are repeatedly variable across datasets for integration
log_print(paste(Sys.time(), 'SelectIntegrationFeatures...'))
spleen_features <- SelectIntegrationFeatures(
    object.list = c(spleen.naive.seurat, spleen.toxoInfected.seurat)
)
# Scaling features for provided objects
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  

# takes 3 minutes
log_print(paste(Sys.time(), 'FindIntegrationAnchors...'))
spleen_anchors <- FindIntegrationAnchors(
    object.list = c(spleen.naive.seurat, spleen.toxoInfected.seurat),
    anchor.features = spleen_features
)
# Finding all pairwise anchors
#   |                                                  | 0 % ~calculating
# Running CCA
# Warning in UseMethod("depth") :
#   no applicable method for 'depth' applied to an object of class "NULL"
# Warning in UseMethod("depth") :
#   no applicable method for 'depth' applied to an object of class "NULL"
# Merging objects
# Finding neighborhoods
# Finding anchors
#     Found 10577 anchors
# Filtering anchors
#     Retained 3285 anchors
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=02m 48s
# Warning message:
# In CheckDuplicateCellNames(object.list = object.list) :
#   Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.


log_print(paste(Sys.time(), 'Integrating data...'))
spleen_integrated <- IntegrateData(anchorset = spleen_anchors)
# Merging dataset 2 into 1
# Extracting anchors for merged samples
# Finding integration vectors
# Finding integration vector weights
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Integrating data

# NOTE: if you look at your seurat object, the default assay as changed from 'RNA' to 'integrated'
# this can be change anytime using the line below
# this would be the same way you would change between scRNA-seq and scATAC-seq
# DefaultAssay(spleen_integrated) <- "RNA"

# Run the standard workflow for visualization and clustering
spleen_integrated <- ScaleData(spleen_integrated, verbose = FALSE)
spleen_integrated <- RunPCA(spleen_integrated, npcs = 30, verbose = FALSE)
spleen_integrated <- RunUMAP(spleen_integrated, reduction = "pca", dims = 1:30)
# Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
# This message will be shown once per session
# 19:00:55 UMAP embedding parameters a = 0.9922 b = 1.112
# 19:00:55 Read 14930 rows and found 30 numeric columns
# 19:00:55 Using Annoy for neighbor search, n_neighbors = 30
# 19:00:55 Building Annoy index with metric = cosine, n_trees = 50
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# 19:00:56 Writing NN index file to temp file /var/folders/vv/f4ftf3qs34z8f896ktyh7p_00000gq/T//RtmpYU4ldC/file806873eeb2c4
# 19:00:56 Searching Annoy index using 1 thread, search_k = 3000
# 19:00:58 Annoy recall = 100%
# 19:00:59 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
# 19:00:59 Initializing from normalized Laplacian + noise (using irlba)
# 19:01:00 Commencing optimization for 200 epochs, with 649038 positive edges
# Using method 'umap'
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# 19:01:06 Optimization finished

spleen_integrated <- FindNeighbors(spleen_integrated, reduction = "pca", dims = 1:30)
# Computing nearest neighbor graph
# Computing SNN

spleen_integrated <- FindClusters(spleen_integrated, resolution = 0.5)
# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

# Number of nodes: 14930
# Number of edges: 550450

# Running Louvain algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Maximum modularity in 10 random starts: 0.9338
# Number of communities: 21
# Elapsed time: 1 seconds


# ----------------------------------------------------------------------
# Visualizing

log_print(paste(Sys.time(), 'Visualizing...'))

DimPlot(spleen_integrated, reduction = "umap", label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'toxoplasma', 'umap_integrated.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# let's see what proportion of our total cells reside in each cluster
# prop.table(table(Idents(spleen_integrated)))

# remember, we have metadata in this integrated seurat object,
# so you can use this to split your UMAP
DimPlot(spleen_integrated, reduction = "umap", 
        split.by = "treatment", # this facets the plot 
        group.by = "seurat_clusters", # labels the cells with values from your group.by variable
        label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'toxoplasma', 'umap_seurat_clusters.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# plot genes of interest on UMAP
FeaturePlot(spleen_integrated, 
            reduction = "umap", 
            features = 'Sdc1',
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE)
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'toxoplasma', 'feature_plot_sdc1.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

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
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'toxoplasma', 'feature_plot_cd4_cd8.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Leveraging cluster identity in your analysis

log_print(paste(Sys.time(), 'Add cluster identity...'))

# now let's rerun our cluster identification using SingleR
spleen_integrated.sce <- as.SingleCellExperiment(spleen_integrated)

# 358 bulk RNA-seq samples of sorted cell populations
MouseRNAseq.data <- MouseRNAseqData(ensembl = FALSE)
predictions <- SingleR(test=spleen_integrated.sce, assay.type.test=1, 
                       ref=MouseRNAseq.data, labels=MouseRNAseq.data$label.main)

#now add back to singleCellExperiment object (or Seurat objects)
spleen_integrated.sce[["SingleR.labels"]] <- predictions$labels
# scater::plotUMAP(spleen_integrated.sce, colour_by = "SingleR.labels")

spleen_integrated2 <- as.Seurat(spleen_integrated.sce, counts = NULL)
DimPlot(
    spleen_integrated2, reduction = "UMAP", 
    split.by = "treatment", # this facets the plot 
    group.by = "SingleR.labels", # labels the cells with values from your group.by variable
    label = TRUE
)
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'toxoplasma', 'umap_singler_labels.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# If we repeat the steps above for different cell type markers, we get a sense for the following cluster IDs
new.cluster.ids <- c("B cells", "RBCs", "CD8+ T cells", "B cells", "RBCs", "CD4+ T cells", "CD4+ T cells", "Monocytes/Macrophages", "Granulocytes", "Monocytes/Macrophages", "B cells", "Plasma cells", "Monocytes/Macrophages", "Monocytes/Macrophages", "Granulocytes", "CD8+ T cells", "CD8+ T cells", "17", "18", "19", "20") 
names(new.cluster.ids) <- levels(spleen_integrated)
spleen_integrated <- RenameIdents(spleen_integrated, new.cluster.ids)
DimPlot(spleen_integrated, reduction = "umap", 
        split.by = "treatment", # this facets the plot 
        label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'toxoplasma', 'umap_treatment.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# take a look at what you've done
# Idents(spleen_integrated)


# ----------------------------------------------------------------------
# subset seurat object to focus on CD4 T cells

spleen_integrated.CD4.Tcells <- subset(spleen_integrated, idents = "CD4+ T cells")

# you could re-cluster if you want...depends on what you're trying to accomplish
#spleen_integrated.CD4.Tcells <- FindNeighbors(spleen_integrated.CD4.Tcells, dims = 1:10, k.param = 5)
#spleen_integrated.CD4.Tcells <- FindClusters(spleen_integrated.CD4.Tcells)
DimPlot(spleen_integrated.CD4.Tcells, reduction = "umap", label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'toxoplasma', 'umap_cd4_t_cells.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# now we need to switch out 'Idents' to be treatment, rather than cluster
Idents(spleen_integrated.CD4.Tcells) <- spleen_integrated.CD4.Tcells$treatment
inf.vs.naive.markers <- FindMarkers(object = spleen_integrated.CD4.Tcells, 
                                    ident.1 = "infected", 
                                    ident.2 = "naive", 
                                    min.pct = 0)
# Note: If you hash out "Idents", you'll get the following cryptic error message:
# Error in WhichCells.Seurat(object = object, idents = ident.1) :
#   Cannot find the following identities in the object: infected"

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
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'toxoplasma', 'feature_plot_ccl5.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
