## This was Step3_multivariate.R

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
source(file.path(wd, 'R', 'utils.R'))
# library(tidyverse) # too broad
library('dplyr')
library('gt') # publication quality tables
library('DT') # for making interactive tables
library('plotly')
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
# Read in data


# From eda.R
log2.cpm.filtered.norm.df <- readr::read_csv(
    file.path(wd, 'data', 'schistosoma', "filtered_normalized_cpm.csv")
)
columns <- colnames(log2.cpm.filtered.norm.df)
sampleLabels <- columns[c(2:length(columns))]


# Read in study design file
# do we really need this?
study_design <- readr::read_tsv(file.path(wd, 'data', 'schistosoma', "studyDesign.txt"))
group <- study_design$sex
group <- factor(group)
timpoint <- study_design$timpoint


# ----------------------------------------------------------------------
# Hierarchical Clustering

# other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
distance <- dist(t(log2.cpm.filtered.norm), method = "maximum")
# other agg methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
clusters <- hclust(distance, method = "average")

# plot and save dendrogram
if (save==TRUE) {
    # see: http://www.sthda.com/english/wiki/creating-and-saving-graphs-r-base-graphs
    png(file.path(wd, 'figures', 'schistosoma', 'pca', 'dendrogram.png'))
    plot(clusters, labels=sampleLabels)
    dev.off()
}



# ----------------------------------------------------------------------
# Perform PCA

# pca_result
pca_result <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)

# Inspect
# ls(pca_result)
# "center"
# "rotation": how much each gene influenced each PC (aka 'scores')
# "scale"
# "sdev"
# "x": how much each sample influenced each PC (aka 'loadings'), note that these have a magnitude and a direction
# summary(pca_result)

# A screeplot is a standard way to view eigenvalues for each PCA
# plot and save screeplot, basically a histogram
if (save==TRUE) {
    png(file.path(wd, 'figures', 'schistosoma', 'pca', 'screeplot.png'))
    screeplot(pca.res)
    dev.off()
}

# compute percentage variance explained by each principal component
variance <- pca_result$sdev^2  # sdev^2 captures these eigenvalues from the PCA result
pct_variance <- round(variance/sum(variance)*100, 1)


# Plot the first two PCs
pca_result.df <- as_tibble(pca_result$x)
fig <- ggplot(pca_result.df) +
    aes(x=PC1, y=PC2, label=sampleLabels) +
    geom_point(size=4) +
    # geom_label() +
    # stat_ellipse() +
    xlab(paste0("PC1 (",pct_variance[1],"%",")")) + 
    ylab(paste0("PC2 (",pct_variance[2],"%",")")) +
    labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
    # coord_fixed() +
    theme_bw()
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'schistosoma', 'pca', 'pca_plot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Create a PCA 'small multiples' chart

pca_result.df <- pca_result$x[,1:4] %>%
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group)
  
pca.pivot <- pivot_longer(pca_result.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot) +
    aes(x=sample, y=loadings, fill=group) +
    geom_bar(stat="identity") +
    facet_wrap(~PC) +
    labs(title="PCA 'small multiples' plot",
         caption=paste0("produced on ", Sys.time())) +
    theme_bw() +
    coord_flip()
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'schistosoma', 'pca_small_multiples_plot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Examples using dplyr

mydata.df <- log2.cpm.filtered.norm.df %>% 
    mutate(mctl.AVG = (MCtl_LEPZQ_1 + MCtl_LEPZQ_2 + MCtl_LEPZQ_3)/3,
           fctl.AVG = (FCtl_LEPZQ_1 + FCtl_LEPZQ_2 + FCtl_LEPZQ_3)/3,
           LogFC = (mctl.AVG - fctl.AVG)) %>% 
    mutate_if(is.numeric, round, 2)

mydata.sort <- mydata.df %>%
    arrange(desc(LogFC)) %>% 
    select(geneID, LogFC)

# Pick out genes of interest 
mydata.filter <- mydata.df %>%
    filter(
        geneID=="MMP1" | geneID=="GZMB" | geneID=="IL1B" | geneID=="GNLY" | geneID=="IFNG" |
        geneID=="CCL4" | geneID=="PRF1" | geneID=="APOBEC3A" | geneID=="UNC13A"
    ) %>%
    select(geneID, mctl.AVG, fctl.AVG, LogFC) %>%
    arrange(desc(LogFC))

# you can also filter based on any regular expression
mydata.grep <- mydata.df %>%
    filter(grepl('CXCL|IFI', geneID)) %>%
    select(geneID, mctl.AVG, fctl.AVG, LogFC) %>%
    arrange(desc(geneID))


# ----------------------------------------------------------------------
# gt table example

gt(mydata.filter)
mydata.filter %>%
    gt() %>%
    fmt_number(columns=2:4, decimals = 1) %>%
    tab_header(title = md("**Regulators of skin pathogenesis**"),
               subtitle = md("*during cutaneous leishmaniasis*")) %>%
    tab_footnote(
        footnote = "Deletion or blockaid ameliorates disease in mice",
        locations = cells_body(
        columns = geneID,
        rows = c(6, 7))) %>% 
    tab_footnote(
        footnote = "Associated with treatment failure in multiple studies",
        locations = cells_body(
        columns = geneID,
        rows = c(2:9))) %>%
    tab_footnote(
        footnote = "Implicated in parasite control",
        locations = cells_body(
        columns = geneID,
        rows = c(2))) %>%
    tab_source_note(
        source_note = md("Reference: Amorim *et al*., (2019). DOI: 10.1126/scitranslmed.aar3619"))


# ----------------------------------------------------------------------
# Searchable table example

DT::datatable(mydata.df[,c(1,12:14)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         #dom = "Blfrtip", 
                         #buttons = c("copy", "csv", "excel"),
                         lengthMenu = c("10", "25", "50", "100")))

# ----------------------------------------------------------------------
# Plotly example

fig <- ggplot(mydata.df) +
  aes(x=mctl.AVG, y=fctl.AVG, 
      text = paste("Symbol:", geneID)) +
  geom_point(shape=16, size=1) +
  ggtitle("disease vs. healthy") +
  theme_bw()

ggplotly(fig)
