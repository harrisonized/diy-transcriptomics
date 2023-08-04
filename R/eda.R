# EDA plots
# This combines 

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
source(file.path(wd, 'R', 'utils.R'))
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


# args
option_list = list(
    make_option(c("-s", "--save"), default=TRUE, action="store_false", metavar="TRUE",
                type="logical", help="disable if you're troubleshooting and don't want to overwrite your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# ----------------------------------------------------------------------
# Pre-script settings

save = opt['save'][[1]]  # save=FALSE

# Start Log
start_time = Sys.time()
log <- log_open(paste("eda ", start_time, '.log', sep=''))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Concatenate data
# Adapted from Step1_TxImport.R

log_print(paste(Sys.time(), 'Reading files...'))

# get human annotations
Tx <- as_tibble(transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name")))
Tx <- dplyr::rename(Tx, target_id = tx_id)  #need to change first column name to 'target_id'
Tx <- dplyr::select(Tx, "target_id", "gene_name")  #transcript ID needs to be the first column in the dataframe

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
sample_ids = unlist(lapply(filepaths, function(x) basename(dirname(x))))
# targets <- readr::read_tsv(file.path(wd, 'data', 'schistosoma', "studyDesign.txt"))  # don't need this

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
    abundance_files,
    type = "kallisto", 
    tx2gene = Tx, 
    txOut = TRUE,  # doesn't work when this is FALSE
    countsFromAbundance = "lengthScaledTPM",
    ignoreTxVersion = TRUE
)
# update column names, they're not included
colnames(txi$counts) <- unlist(lapply(sample_ids, function(x) paste(x, '_count', sep='')))
colnames(txi$abundance) <- unlist(lapply(sample_ids, function(x) paste(x, '_abundance', sep='')))
colnames(txi$length) <- unlist(lapply(sample_ids, function(x) paste(x, '_length', sep='')))

# Inspect
# head(txi$counts)

# Save
if (save==TRUE) {
    log_print(paste(Sys.time(), 'Saving txi files...'))
    if (!dir.exists(file.path(wd, 'data', 'schistosoma', 'txi'))) {
        dir.create(file.path(wd, 'data', 'schistosoma', 'txi'), recursive=TRUE)
    }
    write.table(txi$counts,
                file.path(wd, 'data', 'schistosoma', 'txi', 'txi_counts.csv'),
                quote=FALSE, col.names=TRUE, row.names=TRUE, sep=',')
    write.table(txi$abundance,
                file.path(wd, 'data', 'schistosoma', 'txi', 'txi_abundances.csv'),
                quote=FALSE, col.names=TRUE, row.names=TRUE, sep=',')
    # write.table(txi$length,
    #             file.path(wd, 'data', 'schistosoma', 'txi', 'txi_lengths.csv'),
    #             quote=FALSE, col.names=TRUE, row.names=TRUE, sep=',')
}


# ----------------------------------------------------------------------
# Plot unfiltered counts

log_print(paste(Sys.time(), 'Plotting unfiltered counts...'))

# Generate summary stats
txi$abundance <- transform(
    txi$abundance,
    SD=rowSds(txi$abundance),
    AVG=rowMeans(txi$abundance),
    MED=rowMedians(txi$abundance)
)

# inspect
# head(txi$abundance)
# colSums(txi$abundance)
# colSums(txi$counts)

# Plot tpm unfiltered, non-normalized data
# Let's expand on the plot above a bit more and take a look at each 'layer' of the ggplot code
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
    ggsave(file.path(wd, 'figures', 'schistosoma', 'fig.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Plot unfiltered counts per million

log_print(paste(Sys.time(), 'Plotting unfiltered, un-normalized cpm...'))

# Create DGElist
myDGEList <- edgeR::DGEList(txi$counts)

# note: this doesn't work
if (FALSE) {
    save(myDGEList, file=file.path(wd, 'data', 'schistosoma', "myDGEList.txt"))
    # load(file=file.path(wd, 'data', 'schistosoma', "myDGEList.txt"))  # example load
}

# Compute counts per million
cpm <- edgeR::cpm(myDGEList) 
# colSums(cpm)
log2.cpm <- cpm(myDGEList, log=TRUE)

# reshape data
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", targets$sample)
log2.cpm.df.pivot <- pivot_longer(
    log2.cpm.df, # dataframe to be pivoted
    cols = colnames(log2.cpm.df[, c(2 :length(log2.cpm.df))]), # column names to be stored as a SINGLE variable
    names_to = "samples", # name of that new variable (column)
    values_to = "expression"
)
# log2.cpm.df.pivot  # check data

# note it is easy to plot this pivoted data
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
    ggsave(file.path(wd, 'figures', 'schistosoma', 'p1.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Plot filtered, non-normalized counts per million

log_print(paste(Sys.time(), 'Plotting filtered, non-normalized cpm...'))

# Filter data
# table(rowSums(myDGEList$counts==0)==10)
keepers <- rowSums(cpm>1)>=5  # filter
myDGEList.filtered <- myDGEList[keepers,]
# dim(myDGEList.filtered)

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", targets$sample)
# pivot this FILTERED data, just as you did earlier
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = colnames(log2.cpm.df[, c(2 :length(log2.cpm.df))]), # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

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
    ggsave(file.path(wd, 'figures', 'schistosoma', 'p2.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Plot filtered, normalized counts per million

log_print(paste(Sys.time(), 'Plotting filtered, normalized cpm...'))

# Normalize data
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", targets$sample)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = colnames(log2.cpm.df[, c(2 :length(log2.cpm.df))]), # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

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
    ggsave(file.path(wd, 'figures', 'schistosoma', 'p3.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# ----------------------------------------------------------------------
# Plot all violins together

log_print(paste(Sys.time(), 'Plotting cowplot, normalized cpm...'))

cowplot::plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'schistosoma', 'cowplot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
