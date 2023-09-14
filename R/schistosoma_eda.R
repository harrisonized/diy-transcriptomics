## Adapted from: Step1_TxImport.R, Step2_dataWrangling
## Concatenates the abundance files in the schistosoma dataset and creates some violin plots
## Note that this contains code duplication. However, that is acceptable here, because this is only EDA

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
# library(tidyverse)  # too broad
library('readr')  # read_tsv
suppressMessages(library('dplyr'))
suppressMessages(library('tidyr'))  # pivot_longer
library('tibble')  # as_tibble
suppressMessages(library('matrixStats'))  # rowSums
library('ggplot2')  # ggplot
library('cowplot')  # plot_grid
library('tximport')  # tx_import
suppressMessages(library('EnsDb.Hsapiens.v86'))
suppressMessages(library('GenomicFeatures'))
suppressMessages(library('edgeR'))  # DGElist, cpm
library('optparse')
library('logr')
source(file.path(wd, 'R', 'utils.R'))  # list_files, filter_list_for_match


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default="data/schistosoma/mapped_reads",
                metavar="data/schistosoma/mapped_reads", type="character",
                help="set the input directory"),

    make_option(c("-f", "--filename"), default="abundance.h5",
                metavar="abundance.h5", type="character",
                help="Choose: 'abundance.h5' or 'abundance.tsv'. '.h5' is faster for importing data"),

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
log <- log_open(paste0("schistosoma_eda-", strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading files...'))

# define data source
# schistosoma dataset has 144 files
files = filter_list_for_match(
    list_files(file.path(wd, opt['input-dir'][[1]])),
    opt['filename'][[1]]  # file_ext
)

# get human annotations
tx2gene_obj <- as_tibble(
    GenomicFeatures::transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
)  
tx2gene_obj <- dplyr::rename(tx2gene_obj, target_id = tx_id)  # change first column name to 'target_id'
tx2gene_obj <- dplyr::select(tx2gene_obj, "target_id", "gene_name")  # transcript ID needs to be the first column in the dataframe

# Read abundance files into a txi object
txi <- tximport(
    files,
    type = "kallisto", 
    tx2gene = tx2gene_obj, 
    txOut = TRUE,  # doesn't work when this is FALSE
    countsFromAbundance = "lengthScaledTPM",
    ignoreTxVersion = TRUE
)

# updates columns to `filename_suffix`, eg. F12h_LE_1_abundance
suffix_for_col = c(
    'counts'='_count',
    'abundance'='_abundance',
    'length'='_length'
)
sample_ids = unlist(lapply(files, function(x) basename(dirname(x))))
for (col in names(suffix_for_col)){
    suffix = suffix_for_col[col][[1]]
    colnames(txi[[col]]) <- unlist(lapply(sample_ids, function(x) paste(x, suffix, sep='')))
}

# save
if (!troubleshooting) {

    log_print(paste(Sys.time(), 'Saving txi files...'))
    
    if (!dir.exists(file.path(wd, dirname(opt['input-dir'][[1]]), 'txi'))) {
        dir.create(file.path(wd, dirname(opt['input-dir'][[1]]), 'txi'), recursive=TRUE)
    }
    write.table(txi[['counts']],
                file.path(wd, dirname(opt['input-dir'][[1]]), 'txi', 'txi_counts.csv'),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')
    write.table(txi[['abundance']],
                file.path(wd, dirname(opt['input-dir'][[1]]), 'txi', 'txi_abundances.csv'),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')

    # don't need this
    # write.table(txi$length,
    #             file.path(wd, dirname(opt['input-dir'][[1]]), 'txi', 'txi_lengths.csv'),
    #             quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')
}


# ----------------------------------------------------------------------
# Unfiltered, non-normalized abundance (tpm)

log_print(paste(Sys.time(), 'Plotting unfiltered, non-normalized abundance...'))

# Add SD, AVG, and MED as columns to the abundances matrix
abundance_matrix = as.matrix(txi[['abundance']][, !(colnames(txi[['abundance']]) %in% c('SD', 'AVG', 'MED'))])
txi[['abundance']] <- transform(
    abundance_matrix,
    SD=rowSds(abundance_matrix),
    AVG=rowMeans(abundance_matrix),
    MED=rowMedians(abundance_matrix)
)

fig <- ggplot(txi[['abundance']]) + 
    aes(x = SD, y = MED) +
    geom_point(shape=16, size=2) +
    geom_smooth(method=lm) +
    geom_hex(show.legend = FALSE) +
    labs(y="Median", x = "Standard deviation",
         title="Transcripts per million (TPM)",
         subtitle="unfiltered, non-normalized data",
         caption="DIYtranscriptomics - Spring 2020") +
    theme_classic() +
    theme_dark() + 
    theme_bw()

if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'eda', 'tpm_unfiltered_nonnormalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Create DGElist from counts

# See: https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/DGEList-class
dge_list <- edgeR::DGEList(txi[['counts']])
if (FALSE) {
    # note: this doesn't work for some reason
    save(dge_list, file=file.path(wd, dirname(opt['input-dir'][[1]]), "dge_list.txt"))
    load(file=file.path(wd, dirname(opt['input-dir'][[1]]), "dge_list.txt"))  # example load
}


# ----------------------------------------------------------------------
# Unfiltered, non-normalized log2 counts (cpm)

log_print(paste(Sys.time(), 'Plotting unfiltered, non-normalized log2 counts...'))

cpm_wide <- as_tibble(cpm(dge_list, log=TRUE), rownames = "geneID")
colnames(cpm_wide) <- c("geneID", sample_ids)

# reshape for plotting
cpm_long <- pivot_longer(
    cpm_wide,
    cols = colnames(cpm_wide[, c(2 :length(cpm_wide))]),
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
       # caption=paste0("produced on ", Sys.time())
       ) +
  theme_bw()

if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'eda', 'cpm_unfiltered_nonnormalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Filtered, non-normalized log2 counts (cpm)

log_print(paste(Sys.time(), 'Plotting filtered, non-normalized log2 counts...'))

dge_subset <- dge_list[rowSums(cpm_wide>1)>=5,]
cpm_wide <- as_tibble(cpm(dge_subset, log=TRUE), rownames = "geneID")
colnames(cpm_wide) <- c("geneID", sample_ids)

# reshape for plotting
cpm_long <- pivot_longer(
    cpm_wide,
    cols = colnames(cpm_wide[, c(2 :length(cpm_wide))]),
    names_to = "samples",
    values_to = "expression"
)

p2 <- ggplot(cpm_long) +
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
        subtitle="filtered, non-normalized",
        caption=paste0("produced on ", Sys.time())) +
    theme_bw()
    # coord_flip()

if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'eda', 'tpm_filtered_nonnormalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Filtered, normalized log2 counts (cpm)

log_print(paste(Sys.time(), 'Plotting filtered, normalized log2 counts...'))

dge_subset_norm <- calcNormFactors(dge_subset, method = "TMM")
cpm_norm_wide <- as_tibble(cpm(dge_subset_norm, log=TRUE), rownames = "gene_ID")
colnames(cpm_norm_wide) <- c("gene_ID", sample_ids)
colnames(cpm_norm_wide) <- c("geneID", sample_ids)
if (!troubleshooting) {
    write.table(
        cpm_norm_wide,
        file.path(wd, dirname(opt['input-dir'][[1]]), 'filtered_normalized_cpm.csv'),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep=','
    )
}

# reshape for plotting
cpm_norm_long <- pivot_longer(
    cpm_norm_wide,
    cols = colnames(cpm_norm_wide[, c(2 :length(cpm_norm_wide))]),
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
    ggsave(file.path(wd, opt['output-dir'][[1]], 'eda', 'tpm_filtered_normalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
#  Plot all violins together

log_print(paste(Sys.time(), 'Plotting cowplot...'))

cowplot::plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)

if (!troubleshooting) {
    ggsave(file.path(wd, opt['output-dir'][[1]], 'eda', 'cowplot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
