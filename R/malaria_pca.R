## Adapted from: Step3_multivariate.R
## Requires filtered_normalized_cpm.csv
## Performs PCA
## This is a gutted version of schistosoma_pca.R

wd = dirname(this.path::here())  # wd = '~/github/R/diy-transcriptomics'
# library(tidyverse)  # too broad
library('readr')
library('tibble')
library('tidyr')
suppressMessages(library('dplyr'))
suppressMessages(library('plotly'))
library('optparse')
library('logr')


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default="data/malaria/filtered_normalized_cpm.csv",
                metavar="data/malaria/filtered_normalized_cpm.csv", type="character",
                help="path/to/filtered_normalized_cpm.csv"),

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
cpm_data <- read_csv(file.path(wd, opt['input-file'][[1]]))
sample_ids <- colnames(cpm_data[, 2:ncol(cpm_data)])

# metadata
study_design <- read_tsv(file.path(wd, opt['metadata'][[1]]))
group <- factor(
    unlist(lapply(study_design['group'], function(x) paste0('timepoint_', x)))
)


# ----------------------------------------------------------------------
# Perform PCA

log_print(paste(Sys.time(), 'Perform PCA clustering...'))

# data
pca_result <- prcomp(
    t(cpm_data[, !names(cpm_data)=='target_id']),
    scale.=FALSE, retx=TRUE
)

# Compute percentage variance explained by each PC
variance <- pca_result[['sdev']]^2
pct_variance <- round(variance/sum(variance)*100, 1) 

# Plot the first two PCs
pca_scores <- as_tibble(pca_result[['x']])
pca_plot <- ggplot(pca_scores) +
    aes(x=PC1, y=PC2, color = group) +
    geom_point(size=4) +
    stat_ellipse() +
    xlab(paste0("PC1 (",pct_variance[1],"%",")")) + 
    ylab(paste0("PC2 (",pct_variance[2],"%",")")) +
    labs(title="PCA plot",
         caption=paste0("produced on ", Sys.time())) +
    coord_fixed() +
    theme_bw()

ggplotly(pca_plot)

if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'pca', 'pca_plot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Create a PCA 'small multiples' chart

log_print(paste(Sys.time(), "Plotting 'small multiples' chart..."))

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

fig <- ggplot(pca_scores_long) +
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


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
