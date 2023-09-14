## Adapted from: malarialab_solution.R
## Concatenates the abundance files in the malaria dataset and creates some violin plots
## This script is almost a duplicate of malaria_eda.R


wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
# library(tidyverse) # too broad
suppressMessages(library('readr'))  # read_tsv
suppressMessages(library('tidyr'))  # pivot_longer
suppressMessages(library('edgeR'))  # DGEList
library('ggplot2')
library('tibble')
library('tximport')
library('optparse')
library('logr')
source(file.path(wd, 'R', 'utils.R'))  # list_files, filter_list_for_match


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default="data/malaria/mapped_reads",
                metavar="data/malaria/mapped_reads", type="character",
                help="set the input directory"),

    make_option(c("-f", "--filename"), default="abundance.h5",
                metavar="abundance.h5", type="character",
                help="Choose: 'abundance.h5' or 'abundance.tsv'. '.h5' is faster for importing data"),

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
log <- log_open(paste0("malaria_eda-", strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading files...'))

# define data source
# malaria dataset has 14 files
files = filter_list_for_match(
    list_files(file.path(wd, opt['input-dir'][[1]])),
    opt['filename'][[1]]  # file_ext
)

# Read abundance files into a txi object
txi_gene <- tximport(
    files,
    type = "kallisto",
    # tx2gene = tx2gene_obj,  # the solution doesn't use this
    txOut = TRUE,  # determines whether your data represented at transcript or gene level
    countsFromAbundance = "lengthScaledTPM",
    ignoreTxVersion = TRUE
)

# update columns
suffix_for_col = c(
    'counts'='_count',
    'abundance'='_abundance',
    'length'='_length'
)
sample_ids = unlist(lapply(files, function(x) basename(dirname(x))))
for (col in names(suffix_for_col)){
    suffix = suffix_for_col[col][[1]]
    colnames(txi_gene[[col]]) <- unlist(lapply(sample_ids, function(x) paste(x, suffix, sep='')))
}

# save
if (!troubleshooting) {
    
    log_print(paste(Sys.time(), 'Saving txi files...'))
    
    if (!dir.exists(file.path(wd, 'data', 'malaria', 'txi'))) {
        dir.create(file.path(wd, 'data', 'malaria', 'txi'), recursive=TRUE)
    }
    write.table(txi_gene$counts,
                file.path(wd, dirname(opt['input-dir'][[1]]), 'txi', 'txi_counts.csv'),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')
    write.table(txi_gene$abundance,
                file.path(wd, dirname(opt['input-dir'][[1]]), 'txi', 'txi_abundances.csv'),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')
    
    # write.table(txi_gene$length,
    #             file.path(wd, dirname(opt['input-dir'][[1]]), 'txi', 'txi_lengths.csv'),
    #             quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')
}

# metadata
targets <- read_tsv(file.path(wd, 'data', 'malaria', "studyDesign.txt"))
sample_ids <- targets$sample


# ----------------------------------------------------------------------
# Create DGElist from counts

dge_list <- DGEList(txi_gene[['counts']])


# ----------------------------------------------------------------------
# Raw log2 counts (cpm)

log_print(paste(Sys.time(), 'Plotting raw log2 counts...'))

cpm_wide <- as_tibble(cpm(dge_list, log=TRUE), rownames = "target_id")
colnames(cpm_wide) <- c("target_id", sample_ids)

# reshape for plotting
cpm_long <- pivot_longer(
    cpm_wide,
    cols = -1, # column names to be stored as a SINGLE variable
    names_to = "samples",
    values_to = "expression"
)

p1 <- ggplot(cpm_long) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'eda', 'cpm_unfiltered_nonnormalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Filtered log2 counts (cpm)

log_print(paste(Sys.time(), 'Plotting filtered log2 counts...'))

dge_subset <- dge_list[rowSums(cpm_wide>1)>=2,]
cpm_wide <- as_tibble(cpm(dge_subset, log=TRUE), rownames = "target_id")
colnames(cpm_wide) <- c("target_id", sample_ids)

# reshape for plotting
cpm_long <- pivot_longer(
    cpm_wide,
    cols = -1,
    names_to = "samples",
    values_to = "expression"
)

p2 <- ggplot(cpm_long) +
    aes(x=samples, y=expression, fill=samples) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    stat_summary(
        fun = "median", 
        geom = "point", 
        shape = 95, 
        size = 10, 
        color = "black", 
        show.legend = FALSE) +
    labs(y="log2 expression", x = "sample",
        title="Log2 Counts per Million (CPM)",
        subtitle="filtered, non-normalized",
        caption=paste0("produced on ", Sys.time())) +
    theme_bw()

if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'malaria', 'eda', 'cpm_filtered_nonnormalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Filtered, normalized log2 counts (cpm)

log_print(paste(Sys.time(), 'Plotting filtered, normalized log2 counts...'))

dge_subset_norm <- calcNormFactors(dge_subset, method = "TMM")
cpm_norm_wide <- as_tibble(cpm(dge_subset_norm, log=TRUE), rownames = "target_id")
colnames(cpm_norm_wide) <- c("target_id", sample_ids)

# required for malaria_pca.R
if (!troubleshooting) {
    log_print(paste(Sys.time(), 'Saving filtered, normalized cpm...'))
    write.table(
        cpm_norm_wide,
        file.path(wd, dirname(opt['input-dir'][[1]]), 'filtered_normalized_cpm.csv'),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep=','
    )
}

# reshape for plotting
cpm_norm_long <- pivot_longer(
    cpm_norm_wide,
    cols = -1,
    names_to = "samples",
    values_to = "expression"
)

p3 <- ggplot(cpm_norm_long) +
    aes(x=samples, y=expression, fill=samples) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    stat_summary(fun = "median", 
                 geom = "point", 
                 shape = 95, 
                 size = 10, 
                 color = "black", 
                 show.legend = FALSE) +
    labs(y="log2 expression", x = "sample",
         title="Log2 Counts per Million (CPM)",
         subtitle="filtered, TMM normalized",
         caption=paste0("produced on ", Sys.time())) +
    theme_bw()

if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'malaria', 'eda', 'tpm_filtered_normalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
#  Plot all violins together

log_print(paste(Sys.time(), 'Plotting cowplot...'))

cowplot::plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'malaria', 'eda', 'cowplot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
