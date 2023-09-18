## Adapted from: Step4_publicData.R, Step5_diff_genes.R

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
# library(tidyverse)  # too broad
library('readr')
library('tibble')
library('tidyr')
library('limma')
library('edgeR')  # DGEList
suppressMessages(library('dplyr'))
library('gt')  # publication quality tables
library('DT')  # for making interactive tables
suppressMessages(library('gplots'))
suppressMessages(library('plotly'))
library('RColorBrewer')
library('optparse')
library('logr')


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default="data/malaria/txi/txi_abundances.csv",
                metavar="data/malaria/txi/txi_abundances.csv", type="character",
                help="path/to/txi_abundances.csv"),

    make_option(c("-m", "--metadata"), default="data/malaria/studyDesign.txt",
                metavar="data/malaria/studyDesign.txt", type="character",
                help="path/to/study_design.txt file"),

    make_option(c("-o", "--output-dir"), default="figures/malaria",
                metavar="figures/malaria", type="character",
                help="set the output directory for the figures"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt['troubleshooting'][[1]]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("malaria_pca-", strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read in data

log_print(paste(Sys.time(), 'Reading data...'))

# From eda.R
txi_counts <- read_csv(file.path(wd, opt['input-file'][[1]]))
sample_ids <- colnames(txi_counts)

# metadata
study_design <- read_tsv(file.path(wd, opt['metadata'][[1]]))
group <- factor(
    unlist(lapply(study_design['group'], function(x) paste0('timepoint_', x)))
)


# ----------------------------------------------------------------------
# Differential gene expression analysis

log_print(paste(Sys.time(), "Performing DEG analysis..."))

dge_list <- DGEList(txi_counts)
dge_subset <- dge_list[rowSums(cpm(dge_list)>1)>=2,]
dge_subset_norm <- calcNormFactors(dge_subset, method = "TMM")


# ----------------------------------------------------------------------
# Fit Bayesian model

log_print(paste(Sys.time(), "Fitting Bayesian model..."))

# Use VOOM function from Limma package to model the mean-variance relationship
design_matrix <- model.matrix(~0+group)
colnames(design_matrix) <- levels(group)

png(file.path(wd, opt['output-dir'][[1]], 'pca', 'voom.png'))
log_cpm <- voom(dge_subset_norm, design_matrix, plot = TRUE)  # need an output
dev.off()

fit <- lmFit(log_cpm, design_matrix)  # fit a linear model to your data

# hardcoding is bad practice
contrast.matrix <- makeContrasts(
    DEGs_08hr = timepoint_8hrs - timepoint_0hrs,
    DEGs_16hr = timepoint_16hrs - timepoint_0hrs,
    DEGs_24hr = timepoint_24hrs - timepoint_0hrs,
    DEGs_32hr = timepoint_32hrs - timepoint_0hrs,
    DEGs_40hr = timepoint_40hrs - timepoint_0hrs,
    DEGs_48hr = timepoint_48hrs - timepoint_0hrs,
    levels=design_matrix
)

fits <- contrasts.fit(fit, contrast.matrix)  # extract the linear model
ebFit <- eBayes(fits)  # get bayesian stats
# write.fit(ebFit, file="lmfit_results.txt")

# plot top hits as a table
top_hits <- as_tibble(
    topTable(ebFit, adjust ="BH", coef=6, number=10, sort.by="logFC"),
    rownames = "gene_id"
)
gt_table_1 <- gt(top_hits)
if (!troubleshooting) {
    gtsave(gt_table_1,
           'gt_example_1.png',
           path=file.path(wd, opt['output-dir'][[1]])
    )
}


# ----------------------------------------------------------------------
# Volcano Plots

log_print(paste(Sys.time(), "Making volcano plot..."))

# rerun toptable with high number of genes selected
top_hits <- as_tibble(
    topTable(ebFit, adjust ="BH", coef=1, number=10000, sort.by="logFC"),
    rownames = "gene_id"
)

# now plot
vplot <- ggplot(top_hits) +
    aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", gene_id)) +
    geom_point(size=2) +
    geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
    geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
    geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
    annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
    annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
    labs(title="Volcano plot",
         subtitle = "Malaria erythrocytic cycle",
         caption=paste0("produced on ", Sys.time())) +
    theme_bw()
ggplotly(vplot)

if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'volcano_plot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Venn Diagram

# take a look at what the results of decideTests looks like

log_print(paste(Sys.time(), "Plotting Venn diagram..."))

# decideTests to pull out the DEGs and make Venn Diagram
results <- decideTests(
    ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=1
)

vennDiagram(results[, 1:4], include="up")

if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'pca', 'venn_diagram.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Data Table

log_print(paste(Sys.time(), "Making data table..."))

# retrieve expression data for your DEGs
colnames(log_cpm[['E']]) <- sample_ids

# hardcoding is bad practice
diff_genes <- as_tibble(
    log_cpm[['E']][
        results[,1] !=0 | results[,3] !=0 | results[,4] !=0 |
        results[,5] !=0 | results[,6] !=0,
    ],
    rownames = "gene_id"
)

# create interactive tables to display your DEGs
datatable(
        diff_genes, 
        extensions = c('KeyTable', "FixedHeader"), 
        caption = 'Table 1: DEGs in cutaneous leishmaniasis',
        options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10,
            lengthMenu = c("10", "25", "50", "100"))
    ) %>%
    formatRound(columns=c(2:11), digits=2)

if (!troubleshooting) {
    # NOTE: this .txt file can be directly used for input into other clustering or network analysis tools
    # (e.g., String, Clust (https://github.com/BaselAbujamous/clust, etc.)
    write_tsv(
        diff_genes,
        file.path(wd, dirname(dirname(opt['input-file'][[1]])), "diff_genes.tsv")
    )
}


# ----------------------------------------------------------------------
# Heatmaps and module identification

log_print(paste(Sys.time(), "Plotting heatmap with dendrogram..."))

myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(
    as.dist(1-cor(t(diff_genes[, 2:ncol(diff_genes)]), method="pearson")),
    method="complete"
)
clustColumns <- hclust(
    as.dist(1-cor(diff_genes[, 2:ncol(diff_genes)], method="spearman")),
     method="complete"
)
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

if (!troubleshooting) {
    # see: http://www.sthda.com/english/wiki/creating-and-saving-graphs-r-base-graphs
    png(file.path(wd, opt['output-dir'][[1]], 'heatmap.png'))
    heatmap.2(
        as.matrix(diff_genes[, 2:ncol(diff_genes)]), 
        Rowv=as.dendrogram(clustRows), 
        Colv=as.dendrogram(clustColumns),
        RowSideColors=module.color,
        col=myheatcolors, scale='row', labRow=NA,
        density.info="none", trace="none",  
        cexRow=1, cexCol=1
        # margins=c(8,20)
    )
    dev.off()
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
