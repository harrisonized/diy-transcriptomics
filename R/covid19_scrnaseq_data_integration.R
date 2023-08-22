# This follows DIY_scRNAseq.R.
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
# Integrate multiple scRNA-seq datasets

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
