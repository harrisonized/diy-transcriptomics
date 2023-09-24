## Adapted from: DIY_scRNAseq.R
## Compare two Seurat objects from the spleen of naive mice versus Toxoplasma gondii infected mice
## Runs FindIntegrationAnchors and IntegrateData, then goes through the standard clustering workflow
## Afterward, there are some examples of cluster assignment

wd = dirname(this.path::here())  # wd = '~/github/R/diy-transcriptomics'
source(file.path(wd, 'R', 'functions', 'utils.R'))  # load_rdata
suppressMessages(library('Seurat'))
library('Matrix')
suppressMessages(library('DropletUtils'))
reticulate::use_condaenv('r-reticulate')  # required for cellassign to access tensorflow through python
library('cellassign')  # takes a few seconds
library('scran')
library('SingleR') # automated cell type annotation ('label transfer') using reference data
suppressMessages(library('celldex')) # a large collection of reference expression datasets with curated cell type labels for use with SingleR package
# library('scater')  # plotUMAP, could not install this
library('tibble')
suppressMessages(library('dplyr'))  # %>% operator
library('rjson')
library('R2HTML')
library('readr')
suppressMessages(library('textshape'))
suppressMessages(library('scales'))
library('ggplot2')
suppressMessages(library('DT'))
library('pheatmap')
library('optparse')
library('logr')


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-1"),
                default='data/toxoplasma/spleen.naive.seurat',
                metavar='data/toxoplasma/spleen.naive.seurat',
                type="character",
                help="first seurat object to compare"),

    make_option(c("-j", "--input-2"),
                default='data/toxoplasma/spleen.toxoInfected.seurat',
                metavar='data/toxoplasma/spleen.toxoInfected.seurat',
                type="character",
                help="second seurat object to compare"),

    make_option(c("-l", "--label-1"),  default='naive', metavar='naive',
                type="character", help="label input-1"),

    make_option(c("-k", "--label-2"), default='infected', metavar='infected',
                type="character", help="label input-2"),

    make_option(c("-o", "--output-dir"), default="figures/toxoplasma",
                metavar="figures/toxoplasma", type="character",
                help="set the output directory for the figures"),

    # I'd prefer not to use this
    # make_option(c("-m", "--metadata"), default="data/toxoplasma/studyDesign.txt",
    #             metavar="data/toxoplasma/studyDesign.txt", type="character",
    #             help="path/to/study_design.txt file"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("toxoplasma_data_integration-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading data...'))

seurat_1 <- load_rdata(file.path(wd, opt[['input-1']]))  # spleen.naive.seurat
seurat_1$treatment <- opt[['label-1']]  # naive
# Note: I normally don't like to use the '$' operator, but in this case, 
# seurat_1$treatment adds the label to seurat_1@meta.data[['treatment']], and
# seurat_1[['treatment']] does not

# # inspect
# seurat_1@meta.data[['treatment']]

# Plot UMAP of raw data
DimPlot(seurat_1, reduction = "umap", split.by = "orig.ident", label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], 'umap-unlabeled_1.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

seurat_2 <- load_rdata(file.path(wd, opt[['input-2']]))  # spleen.toxoInfected.seurat
seurat_2$treatment <- opt[['label-2']]  # infected

# Plot UMAP of raw data
DimPlot(seurat_2, reduction = "umap", split.by = "orig.ident", label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], 'umap-unlabeled_2.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Integrate data

log_print(paste(Sys.time(), 'Integrating data...'))

# select features that are repeatedly variable across datasets for integration
integration_features <- SelectIntegrationFeatures(object.list = c(seurat_1, seurat_2))

# takes 3 minutes
integration_anchors <- FindIntegrationAnchors(
    object.list = c(seurat_1, seurat_2),
    anchor.features = integration_features
)

# Scaling features for provided objects
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
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
#   Some cell names are duplicated across objects provided. Renaming to enforce
#   unique cell names.

integrated_seurat <- IntegrateData(anchorset = integration_anchors)
# Merging dataset 2 into 1
# Extracting anchors for merged samples
# Finding integration vectors
# Finding integration vector weights
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Integrating data

# NOTE: in the seurat object, the default assay as changed from 'RNA' to 'integrated'
# this can be change anytime using the line below
# this would be the same way you would change between scRNA-seq and scATAC-seq
# DefaultAssay(integrated_seurat) <- "RNA"


# ----------------------------------------------------------------------
# Standard clustering workflow

log_print(paste(Sys.time(), 'Running standard clustering workflow...'))

integrated_seurat <- ScaleData(integrated_seurat, verbose = FALSE)
integrated_seurat <- RunPCA(integrated_seurat, npcs = 30, verbose = FALSE)
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "pca", dims = 1:30)
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

integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "pca", dims = 1:30)
# Computing nearest neighbor graph
# Computing SNN

integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5)
# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#
# Number of nodes: 14930
# Number of edges: 550450
#
# Running Louvain algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Maximum modularity in 10 random starts: 0.9338
# Number of communities: 21
# Elapsed time: 1 seconds

# check proportion of total cells that reside in each cluster
# prop.table(table(Idents(integrated_seurat)))

# Plot UMAP with clusters highlighted
DimPlot(integrated_seurat, reduction = "umap", label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], 'umap-integrated.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# Split UMAP by treatment
DimPlot(integrated_seurat, reduction = "umap", 
        split.by = "treatment",  # facets the plot 
        group.by = "seurat_clusters", # labels the cells with values from your group.by variable
        label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], 'umap-integrated-split.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# Plot UMAP for specific gene of interest
gene = 'Sdc1'
FeaturePlot(integrated_seurat, 
            reduction = "umap", 
            features = gene,
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], paste0('umap-', tolower(gene), '.png')),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# Plot UMAP for multiple genes of interest
genes <- c("Cd4", "Cd8a")
FeaturePlot(integrated_seurat, 
            reduction = "umap", 
            features = genes,
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']],
                     paste0('umap-', paste(tolower(genes), collapse="_"), '.png')),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Label clusters using celldex

log_print(paste(Sys.time(), 'Labeling clusters...'))

sce_counts <- as.SingleCellExperiment(integrated_seurat)
ref_data <- MouseRNAseqData(ensembl = FALSE)
predictions <- SingleR(
    test=sce_counts,
    assay.type.test=1, 
    ref=ref_data,
    labels=ref_data[['label.main']]
)

# now add back to singleCellExperiment object (or Seurat objects)
sce_counts$SingleR.labels <- predictions[['labels']]
# scater::plotUMAP(sce_counts, colour_by = "SingleR.labels")

seurat_counts <- as.Seurat(sce_counts, counts = NULL)
DimPlot(
    seurat_counts, reduction = "UMAP", 
    split.by = "treatment", # this facets the plot 
    group.by = "SingleR.labels", # labels the cells with values from your group.by variable
    label = TRUE
)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], 'labeled', 'umap-integrated-celldex_labeled.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Manually assign cluster labels

# If we repeat the steps above for different cell type markers, we get a sense for the following cluster IDs
new.cluster.ids <- c(
    "B cells", "RBCs", "CD8+ T cells", "B cells", "RBCs", "CD4+ T cells", "CD4+ T cells",
    "Monocytes/Macrophages", "Granulocytes", "Monocytes/Macrophages",
    "B cells", "Plasma cells", "Monocytes/Macrophages", "Monocytes/Macrophages",
    "Granulocytes", "CD8+ T cells", "CD8+ T cells", "17", "18", "19", "20"
) 
names(new.cluster.ids) <- levels(integrated_seurat)
integrated_seurat <- RenameIdents(integrated_seurat, new.cluster.ids)
DimPlot(integrated_seurat, reduction = "umap", 
        split.by = "treatment", # this facets the plot 
        label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], 'labeled', 'umap-integrated-manually_labeled.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# inspect
# Idents(integrated_seurat)


# ----------------------------------------------------------------------
# Find cluster-specific DEGs

cluster_of_interest = "CD4+ T cells"
integrated_seurat_subset <- subset(integrated_seurat, idents = cluster_of_interest)

# you could re-cluster if you want...depends on what you're trying to accomplish
#integrated_seurat_subset <- FindNeighbors(integrated_seurat_subset, dims = 1:10, k.param = 5)
#integrated_seurat_subset <- FindClusters(integrated_seurat_subset)

DimPlot(integrated_seurat_subset, reduction = "umap", label = TRUE)
if (!troubleshooting) {

    # format string
    cluster_of_interest_f <- unlist(strsplit(cluster_of_interest, split=' ')) %>%
        paste(collapse="_") %>%
        gsub('[+]', '', .) %>%  # remove special characters
        tolower()

    ggsave(file.path(wd, opt[['output-dir']], 'labeled',
                     paste0('umap-', cluster_of_interest_f, '.png')),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# Switch out 'Idents' to be treatment, rather than cluster
# Note: If this line is not run, then you will get the following error when running FindMarkers:
# Error in WhichCells.Seurat(object = object, idents = ident.1) :
#   Cannot find the following identities in the object: infected"
Idents(integrated_seurat_subset) <- integrated_seurat_subset[['treatment']]

degs <- FindMarkers(
    object = integrated_seurat_subset, 
    ident.1 = "infected", 
    ident.2 = "naive", 
    min.pct = 0
)

degs[['pct.diff']] <- degs[['pct.1']] - degs[['pct.2']]
degs.df <- as_tibble(degs, rownames = "geneID")

# # top 20 DEGs
# degs_top20 <- degs.df %>% arrange(desc(avg_log2FC)) %>% slice(1:20)
# if (!troubleshooting) {
#     write_tsv(
#         degs_top20,
#         file.path(wd, opt[['input-dir']], 'degs-top_20.tsv')
#     )
# }

# Plot UMAP for specific gene of interest
gene = "Ccl5"
FeaturePlot(integrated_seurat_subset, 
            reduction = "umap", 
            features = gene,
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], 'labeled',
                     paste0('umap-', tolower(gene), '.png')),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()


# ----------------------------------------------------------------------
# Code graveyard

# # example of annotating seurat objects using metadata
# study_design <- read_tsv(file.path(wd, opt[['metadata']]))
# seurat_1$treatment <- study_design[['treatment']][1]  # naive
# seurat_2$treatment <- study_design[['treatment']][2]  # infected
