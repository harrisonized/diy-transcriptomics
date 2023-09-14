## Source: Step1_TxImport.R, Step2_dataWrangling
## 1. Concat data
## 2. Filter and normalization
## 3. Plot

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
# library(tidyverse)  # too broad
library('readr')  # read_tsv
library('tibble')  # as_tibble
library('tidyr')  # pivot_longer
library('matrixStats')  # colSums, rowSums
library('ggplot2')  # ggplot
library('cowplot')  # plot_grid

library('tximport')  # tx_import
suppressMessages(library('EnsDb.Hsapiens.v86'))
suppressMessages(library('edgeR'))  # DGElist, cpm
library('optparse')
library('logr')
source(file.path(wd, 'R', 'utils.R'))  # list_files, filter_list_for_match


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(

    make_option(c("-i", "--input-dir"), default="data/schistosoma/mapped_reads",
                metavar="data/leishmania/mapped_reads", type="character",
                help="set the input directory"),

    make_option(c("-f", "--file-ext"), default="h5",
                metavar="h5", type="character",
                help="Choose: 'h5' or 'tsv'. 'h5' is faster for importing data"),

    make_option(c("-s", "--save"), default=TRUE, action="store_false", metavar="TRUE",
                type="logical", help="disable if you're troubleshooting and don't want to overwrite your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

save = opt['save'][[1]]  # save=FALSE

# Start Log
start_time = Sys.time()
log <- log_open(paste0("schistosoma_eda-", strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Concat data
# Adapted from Step1_TxImport.R

log_print(paste(Sys.time(), 'Reading files...'))

# get human annotations
tx2gene_obj <- as_tibble(
    GenomicFeatures::transcripts('EnsDb.Hsapiens.v86', columns=c("tx_id", "gene_name"))
)
tx2gene_obj <- dplyr::rename(tx2gene_obj, target_id = tx_id)  #need to change first column name to 'target_id'
tx2gene_obj <- dplyr::select(tx2gene_obj, "target_id", "gene_name")  #transcript ID needs to be the first column in the dataframe

# Example file: "data/schistosoma/mapped_reads/F12h_LE_1/abundance.tsv"
#   target_id    length eff_length est_counts   tpm
#   <chr>         <dbl>      <dbl>      <dbl> <dbl>
# 1 Smp_186980.1    402      239.       171   58.5 
# 2 Smp_197050.1    186       50.3       59   95.8 
# 3 Smp_160500.1   5316     5152.       244    3.87
abundance_filepaths = filter_list_for_match(
    list_files(file.path(wd, "data", "schistosoma", "mapped_reads")),
    'abundance.h5'
)

# Read abundance files into a txi object
# see: https://github.com/mikelove/tximport/blob/devel/R/tximport.R
# names(txi)
# [1] "abundance"  "counts"  "infReps"  "length"  "countsFromAbundance
# txi$abundance = tpm
# txi$counts != counts, they did something
# txi$infReps: inferential replicates, not sure what this is
# txi$length = eff_length
# txi$countsFromAbundance = "lengthScaledTPM"
txi <- tximport(
    abundance_filepaths,
    type = "kallisto", 
    tx2gene = tx2gene_obj, 
    txOut = TRUE,  # doesn't work when this is FALSE
    countsFromAbundance = "lengthScaledTPM",
    ignoreTxVersion = TRUE
)


# update column names, they're not included
sample_ids = unlist(lapply(abundance_filepaths, function(x) basename(dirname(x))))
colnames(txi$counts) <- unlist(lapply(sample_ids, function(x) paste(x, '_count', sep='')))
colnames(txi$abundance) <- unlist(lapply(sample_ids, function(x) paste(x, '_abundance', sep='')))
colnames(txi$length) <- unlist(lapply(sample_ids, function(x) paste(x, '_length', sep='')))

if (save==TRUE) {
    log_print(paste(Sys.time(), 'Saving txi files...'))
    if (!dir.exists(file.path(wd, 'data', 'schistosoma', 'txi'))) {
        dir.create(file.path(wd, 'data', 'schistosoma', 'txi'), recursive=TRUE)
    }
    write.table(txi$counts,
                file.path(wd, 'data', 'schistosoma', 'txi', 'txi_counts.csv'),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')
    write.table(txi$abundance,
                file.path(wd, 'data', 'schistosoma', 'txi', 'txi_abundances.csv'),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')
    # write.table(txi$length,
    #             file.path(wd, 'data', 'schistosoma', 'txi', 'txi_lengths.csv'),
    #             quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')
}

# Add SD, AVG, and MED as columns to the txi$abundance tibble for plotting later
txi$abundance <- transform(
    txi$abundance,
    SD=rowSds(txi$abundance),
    AVG=rowMeans(txi$abundance),
    MED=rowMedians(txi$abundance)
)


# ----------------------------------------------------------------------
# Transform data
# Adapted from Step2_dataWrangling


# Create DGElist
# See: https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/DGEList-class
dge_list <- edgeR::DGEList(txi$counts)

# note: this doesn't work for some reason
if (FALSE) {
    save(dge_list, file=file.path(wd, 'data', 'schistosoma', "dge_list.txt"))
    load(file=file.path(wd, 'data', 'schistosoma', "dge_list.txt"))  # example load
}

# unfiltered, non-normalized cpm
cpm <- edgeR::cpm(dge_list)
log2.cpm <- edgeR::cpm(dge_list, log=TRUE)
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sample_ids)
log2.cpm.df.pivot <- pivot_longer(
    log2.cpm.df, # dataframe to be pivoted
    # column names to be stored as a SINGLE variable
    cols = colnames(log2.cpm.df[, c(2 :length(log2.cpm.df))]),
    names_to = "samples", # name of that new variable (column)
    values_to = "expression"
)


# filtered, non-normalized cpm
# table(rowSums(dge_list$counts==0)==10)
keepers <- rowSums(cpm>1)>=5  # filter
dge_list.filtered <- dge_list[keepers,]
# dim(dge_list.filtered)

# filtered, non-normalized cpm
log2.cpm.filtered <- cpm(dge_list.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sample_ids)
# pivot this FILTERED data, just as you did earlier
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = colnames(log2.cpm.df[, c(2 :length(log2.cpm.df))]), # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data) n

# filtered, normalized cpm
# need this for the next script
dge_list.filtered.norm <- calcNormFactors(dge_list.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(dge_list.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "gene_ID")
colnames(log2.cpm.filtered.norm.df) <- c("gene_ID", sample_ids)

if (save==TRUE) {
    log_print(paste(Sys.time(), 'Saving filtered, normalized cpm...'))
    write.table(
        log2.cpm.filtered.norm.df,
        file.path(wd, 'data', 'schistosoma', 'filtered_normalized_cpm.csv'),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep=','
    )
}


log2.cpm.filtered.norm.df.pivot <- pivot_longer(
    log2.cpm.filtered.norm.df, # dataframe to be pivoted
    cols = colnames(log2.cpm.df[, c(2 :length(log2.cpm.df))]), # column names to be stored as a SINGLE variable
    names_to = "samples", # name of that new variable (column)
    values_to = "expression"
) # name of new variable (column) storing all the values (data)


# ----------------------------------------------------------------------
# Plot
# Adapted from Step2_dataWrangling


# Plot tpm unfiltered, non-normalized tpm
log_print(paste(Sys.time(), 'Plotting unfiltered, non-normalized tpm...'))
fig <- ggplot(txi$abundance) + 
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

if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'schistosoma', 'eda', 'tpm_unfiltered_nonnormalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# Plot unfiltered, un-normalized log2 cpm
log_print(paste(Sys.time(), 'Plotting unfiltered, un-normalized log2 cpm...'))
p1 <- ggplot(log2.cpm.df.pivot) +
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
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'schistosoma', 'eda', 'cpm_unfiltered_nonnormalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# Plot filtered, non-normalized counts per million
log_print(paste(Sys.time(), 'Plotting filtered, non-normalized cpm...'))
p2 <- ggplot(log2.cpm.filtered.df.pivot) +
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
# Also: try using coord_flip() at the end of the ggplot code

if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'schistosoma', 'eda', 'tpm_filtered_nonnormalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# Plot filtered, normalized cpm
log_print(paste(Sys.time(), 'Plotting filtered, normalized cpm...'))
p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
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
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'schistosoma', 'eda', 'tpm_filtered_normalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# Plot all violins together
log_print(paste(Sys.time(), 'Plotting cowplot...'))
cowplot::plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'schistosoma', 'eda', 'cowplot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
