## This was Step3_multivariate.R

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
source(file.path(wd, 'R', 'utils.R'))
# library(tidyverse) # too broad
library('readr')
library('tibble')
library('tidyr')
library('limma')
library('edgeR')  # DGEList
library('dplyr')
library('gt') # publication quality tables
library('DT') # for making interactive tables
library('gplots')
library('plotly')
library('RColorBrewer')
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
log <- log_open(paste("pca ", start_time, '.log', sep=''))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read in data

# From eda.R
log2.cpm.filtered.norm <- readr::read_csv(
    file.path(wd, 'data', 'malaria', "filtered_normalized_cpm.csv"),
)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
columns <- colnames(log2.cpm.filtered.norm)
sampleLabels <- columns[c(1:length(columns))]

txi_counts <- readr::read_csv(
    file.path(wd, 'data', 'malaria', 'txi', 'txi_abundances.csv'),
)

# study_design <- readr::read_tsv(file.path(wd, 'data', 'malaria', "studyDesign.txt"))
# this is bad practice
group <- factor(c("timepoint_00hr",
                  "timepoint_00hr",
                  "timepoint_08hr",
                  "timepoint_08hr",
                  "timepoint_16hr", 
                  "timepoint_16hr",
                  "timepoint_24hr", 
                  "timepoint_24hr",
                  "timepoint_32hr", 
                  "timepoint_32hr",
                  "timepoint_40hr",
                  "timepoint_40hr",
                  "timepoint_48hr",
                  "timepoint_48hr"))


# ----------------------------------------------------------------------
# Perform PCA

# data
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2  # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)

# plot
pca.plot <- ggplot(pca.res.df) +
    aes(x=PC1, y=PC2, color = group) +
    geom_point(size=4) +
    stat_ellipse() +
    xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
    ylab(paste0("PC2 (",pc.per[2],"%",")")) +
    labs(title="PCA plot",
         caption=paste0("produced on ", Sys.time())) +
    coord_fixed() +
    theme_bw()

ggplotly(pca.plot)

if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'malaria', 'pca', 'pca_plot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Create a PCA 'small multiples' chart

# this is another way to view PCA laodings to understand impact of each sample on each pricipal component
pca.res.df <- pca.res$x[,1:4] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
    as_tibble() %>%
    add_column(sample = sampleLabels,
               group = group)

pca.pivot <- pivot_longer(
    pca.res.df,
    cols = PC1:PC4, # column names to be stored as a SINGLE variable
    names_to = "PC", # name of that new variable (column)
    values_to = "loadings"
)  # name of new variable (column) storing all the values (data)

fig <- ggplot(pca.pivot) +
    aes(x=sample, y=loadings, fill=group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
    geom_bar(stat="identity") +
    facet_wrap(~PC) +
    labs(title="PCA 'small multiples' plot",
         caption=paste0("produced on ", Sys.time())) +
    theme_bw() +
    coord_flip()

if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'malaria', 'pca', 'pca_small_multiples_plot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Differential gene expression analysis (step 5 script) ----

# group <- factor(targets$group)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

# Model mean-variance trend and fit linear model to data
# Use VOOM function from Limma package to model the mean-variance relationship
myDGEList <- DGEList(txi_counts)
keepers <- rowSums(cpm(myDGEList)>1)>=2  #user defined
myDGEList.filtered <- myDGEList[keepers,]
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)


# ----------------------------------------------------------------------
# Fit Bayesian model

# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)

# Contrast matrix
contrast.matrix <- makeContrasts(
    DEGs_08hr = timepoint_08hr - timepoint_00hr,
    DEGs_16hr = timepoint_16hr - timepoint_00hr,
    DEGs_24hr = timepoint_24hr - timepoint_00hr,
    DEGs_32hr = timepoint_32hr - timepoint_00hr,
    DEGs_40hr = timepoint_40hr - timepoint_00hr,
    DEGs_48hr = timepoint_48hr - timepoint_00hr,
    levels=design
)

#  extract the linear model fit
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)
#write.fit(ebFit, file="lmfit_results.txt")

# TopTable to view DEGs
myTopHits <- topTable(ebFit, adjust ="BH", coef=6, number=10, sort.by="logFC")

# convert to a tibble
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

gt_table_1 <- gt(myTopHits.df)

gtsave(gt_table_1,
       'gt_example_1.png',
       path=file.path(wd, 'figures', 'malaria')
)

# ----------------------------------------------------------------------
# Volcano Plots

# rerun toptable with high number of genes selected
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=10000, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")
# now plot
vplot <- ggplot(myTopHits.df) +
    aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
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
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'malaria', 'volcano_plot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Venn Diagram

# decideTests to pull out the DEGs and make Venn Diagram
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=1)

# take a look at what the results of decideTests looks like
# head(results)
# summary(results)
vennDiagram(results[, 1:4], include="up")
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'malaria', 'pca', 'venn_diagram.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Data Table

# retrieve expression data for your DEGs
# head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0 | results[,3] !=0 | results[,4] !=0 | results[,5] !=0 | results[,6] !=0,]
# head(diffGenes)
# dim(diffGenes)
#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

# create interactive tables to display your DEGs
datatable(
        diffGenes.df, 
        extensions = c('KeyTable', "FixedHeader"), 
        caption = 'Table 1: DEGs in cutaneous leishmaniasis',
        options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10,
            lengthMenu = c("10", "25", "50", "100"))
    ) %>%
    formatRound(columns=c(2:11), digits=2)

if (save==TRUE) {
    # NOTE: this .txt file can be directly used for input into other clustering or network analysis tools
    # (e.g., String, Clust (https://github.com/BaselAbujamous/clust, etc.)
    write_tsv(diffGenes.df, file.path(wd, 'data', 'malaria', "DiffGenes.txt"))
}


# ----------------------------------------------------------------------
# Heatmaps and module identification

myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 


# plot and save dendrogram
if (save==TRUE) {
    # see: http://www.sthda.com/english/wiki/creating-and-saving-graphs-r-base-graphs
    png(file.path(wd, 'figures', 'malaria', 'heatmap.png'))
    heatmap.2(
        diffGenes, 
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


# ----------------------------------------------------------------------
# Code graveyard

# omg this is so bad
# mydata.df <- mutate(
#         log2.cpm.filtered.norm.df,
#         timepoint_0hr = (totalRNA_0hr_rep1 + totalRNA_0hr_rep2)/2, 
#         timepoint_8hr = (totalRNA_8hr_rep1 + totalRNA_8hr_rep2)/2, 
#         timepoint_16hr = (totalRNA_16hr_rep1 + totalRNA_16hr_rep2)/2, 
#         timepoint_24hr = (totalRNA_24hr_rep1 + totalRNA_24hr_rep2)/2, 
#         timepoint_32hr = (totalRNA_32hr_rep1 + totalRNA_32hr_rep2)/2, 
#         timepoint_48hr = (totalRNA_48hr_rep1 + totalRNA_48hr_rep2)/2, 
#         #now make columns comparing each of the averages above that you're interested in
#         LogFC = (timepoint_32hr - timepoint_0hr)
#     ) %>% #note that this is the first time you've seen the 'pipe' operator
#     mutate_if(is.numeric, round, 2)

# DT::datatable(
#     mydata.df[,c(1,16,20,22)], 
#     extensions = c('KeyTable', "FixedHeader"), 
#     filter = 'top',
#     options = list(
#         keys = TRUE, 
#         searchHighlight = TRUE, 
#         pageLength = 10, 
#         lengthMenu = c("10", "25", "50", "100")
#     )
# )
