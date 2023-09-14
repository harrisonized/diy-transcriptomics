## Adapted from: Step3_multivariate.R
## Requires filtered_normalized_cpm.csv
## Performs PCA and plots the dendrogram
## Then, there are examples from gt, DT, and plotly

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
# library(tidyverse) # too broad
library('tibble')
library('tidyr')
suppressMessages(library('dplyr'))
library('gt')  # publication quality tables
library('DT')  # interactive tables
suppressMessages(library('plotly'))
library('optparse')
library('logr')


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default="data/schistosoma/filtered_normalized_cpm.csv",
                metavar="data/schistosoma/filtered_normalized_cpm.csv", type="character",
                help="path/to/filtered_normalized_cpm.csv"),

    make_option(c("-o", "--output-dir"), default="figures/schistosoma",
                metavar="figures/schistosoma", type="character",
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
log <- log_open(paste0("schistosoma_pca-", strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read in data

# cpm
cpm_data <- readr::read_csv(file.path(wd, opt['input-file'][[1]]))
sample_ids <- colnames(cpm_data[, 2:ncol(cpm_data)])

study_design <- readr::read_tsv(file.path(wd, 'data', 'schistosoma', "studyDesign.txt"))
group <- factor(study_design[['sex']])


# ----------------------------------------------------------------------
# Hierarchical Clustering

log_print(paste(Sys.time(), 'Hierarchical clustering...'))

# distance methods: "euclidean", maximum", "manhattan", "canberra", "binary", "minkowski"
# agg methods: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
distances <- dist(
    t(cpm_data[, !names(cpm_data)=='target_id']),
    method = "maximum"
)
clusters <- hclust(distances, method = "average")

# plot dendrogram
if (!troubleshooting) {
    # see: http://www.sthda.com/english/wiki/creating-and-saving-graphs-r-base-graphs
    png(file.path(wd, opt['output-dir'][[1]], 'pca', 'dendrogram.png'))
    plot(clusters, labels=sample_ids)
    dev.off()
}


# ----------------------------------------------------------------------
# Perform PCA

log_print(paste(Sys.time(), 'Perform PCA clustering...'))

# pca_result
# See: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp
# rotation: how much each sample influenced each PC (aka 'loadings')
# x: how much each gene influenced each PC (aka 'scores')
# sdev
# center, scale
# summary(pca_result)
pca_result <- prcomp(
    t(cpm_data[, !names(cpm_data)=='target_id']),
    scale.=FALSE, retx=TRUE
)

# Inspect
# ls(pca_result)

# Plot eigenvalues for each PC
if (!troubleshooting) {
    png(file.path(wd, opt['output-dir'][[1]], 'pca', 'screeplot.png'))
    screeplot(pca_result)
    dev.off()
}

# Compute percentage variance explained by each PC
variance <- pca_result[['sdev']]^2  # sdev^2 captures these eigenvalues from the PCA result
pct_variance <- round(variance/sum(variance)*100, 1)

# Plot the first two PCs
pca_scores <- as_tibble(pca_result[['x']])
fig <- ggplot(pca_scores) +
    aes(x=PC1, y=PC2, label=sample_ids) +
    geom_point(size=4) +
    # geom_label() +
    # stat_ellipse() +
    xlab(paste0("PC1 (",pct_variance[1],"%",")")) + 
    ylab(paste0("PC2 (",pct_variance[2],"%",")")) +
    labs(title="PCA plot",
         caption=paste0("produced on ", Sys.time())) +
    # coord_fixed() +
    theme_bw()

if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'pca', 'pca_plot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Plot "small multiples" chart

pca_scores <- pca_result[['x']][,1:4] %>%
  as_tibble() %>%
  add_column(sample = sample_ids,
             group = group)

# reshape for plotting
pca_scores_long <- pivot_longer(
    pca_scores,
    cols = PC1:PC4,
    names_to = "PC",
    values_to = "loadings"
)

fig2 <- ggplot(pca_scores_long) +
    aes(x=sample, y=loadings, fill=group) +
    geom_bar(stat="identity") +
    facet_wrap(~PC) +
    labs(title="PCA 'small multiples' plot",
         caption=paste0("produced on ", Sys.time())) +
    theme_bw() +
    coord_flip()
if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'pca', 'pca_small_multiples_plot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Examples using dplyr

cpm_agg <- cpm_data %>% 
    mutate(male_ctl_avg = (MCtl_LEPZQ_1 + MCtl_LEPZQ_2 + MCtl_LEPZQ_3)/3,
           female_ctl_avg = (FCtl_LEPZQ_1 + FCtl_LEPZQ_2 + FCtl_LEPZQ_3)/3,
           log_fold_change = (male_ctl_avg - female_ctl_avg)) %>% 
    mutate_if(is.numeric, round, 2) %>%
    arrange(desc(log_fold_change)) %>% 
    dplyr::select(target_id, male_ctl_avg, female_ctl_avg, log_fold_change)

# Example of filtering
target_ids <- c("Smp_017610.1", "Smp_197050.1", "Smp_160500.1", "Smp_160490.1", "Smp_070540.1")
cpm_subset <- cpm_agg %>%
    dplyr::filter(target_id %in% target_ids) %>%
    dplyr::select(target_id, male_ctl_avg, female_ctl_avg, log_fold_change) %>%
    arrange(desc(log_fold_change))

# you can also filter based on any regular expression
# cpm_subset <- cpm_agg %>%
#     dplyr::filter(grepl('Smp_017[[:digit:]]+', target_id)) %>%
#     dplyr::select(target_id, male_ctl_avg, female_ctl_avg, log_fold_change) %>%
#     arrange(desc(target_id))


# ----------------------------------------------------------------------
# gt table example

gt_table_1 <- gt(cpm_subset)
gtsave(gt_table_1,
       'gt_example_1.png',
       path=file.path(wd, opt['output-dir'][[1]])
)

gt_table_2 <- cpm_subset %>%
    gt() %>%
    fmt_number(columns=1:3, decimals = 1) %>%
    tab_header(title = md("**Regulators of skin pathogenesis**"),
               subtitle = md("*during cutaneous leishmaniasis*")) %>%
    tab_footnote(
        footnote = "Deletion or blockaid ameliorates disease in mice",
        locations = cells_body(
        columns = target_id,
        rows = c(4, 5))) %>% 
    tab_footnote(
        footnote = "Associated with treatment failure in multiple studies",
        locations = cells_body(
        columns = target_id,
        rows = c(2:5))) %>%
    tab_footnote(
        footnote = "Implicated in parasite control",
        locations = cells_body(
        columns = target_id,
        rows = c(2))) %>%
    tab_source_note(
        source_note = md("Reference: Amorim *et al*., (2019). DOI: 10.1126/scitranslmed.aar3619"))
if (!troubleshooting) {
    gtsave(gt_table_2,
           'gt_example_2.png',
           path=file.path(wd, opt['output-dir'][[1]])
    )
}


# ----------------------------------------------------------------------
# Searchable table example

DT::datatable(
    cpm_agg, 
    extensions = c('KeyTable', "FixedHeader"), 
    filter = 'top',
    options = list(keys = TRUE, 
                   searchHighlight = TRUE, 
                   pageLength = 10, 
                   #dom = "Blfrtip", 
                   #buttons = c("copy", "csv", "excel"),
                   lengthMenu = c("10", "25", "50", "100"))
)
# Warning message:
# In instance$preRenderHook(instance) :
#   It seems your data is too big for client-side DataTables.
#   You may consider server-side processing: https://rstudio.github.io/DT/server.html


# ----------------------------------------------------------------------
# Plotly example

fig <- ggplot(cpm_agg) +
  aes(x=male_ctl_avg, y=female_ctl_avg, 
      text = paste("Symbol:", target_id)) +
  geom_point(shape=16, size=1) +
  ggtitle("disease vs. healthy") +
  theme_bw()

ggplotly(fig)

if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'pca', 'plotly_scatter_example.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
