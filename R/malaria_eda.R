# Lab 9
# Ask is to go through the Step_1 through Step_6 scripts and generate some figures


wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
source(file.path(wd, 'R', 'utils.R'))
# library(tidyverse) # too broad
library('readr')  # read_tsv
library('tidyr')  # pivot_longer
library('edgeR')  # DGEList
library('ggplot2')
library('tibble')
library('tximport')
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
log <- log_open(paste0("malaria_eda-", strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

# metadata
targets <- read_tsv(file.path(wd, 'data', 'malaria', "studyDesign.txt"))
sample_ids <- targets$sample

# main data
paths = filter_list_for_match(
    list_files(file.path(wd, 'data', 'malaria')),
    'abundance.h5'
)
txi_gene <- tximport(
    paths,
    type = "kallisto",
    txOut = TRUE, # determines whether your data represented at transcript or gene level
    countsFromAbundance = "lengthScaledTPM",
    ignoreTxVersion = TRUE
)
# update column names, they're not included
sample_ids = unlist(lapply(paths, function(x) basename(dirname(x))))
colnames(txi_gene$counts) <- unlist(lapply(sample_ids, function(x) paste(x, '_count', sep='')))
colnames(txi_gene$abundance) <- unlist(lapply(sample_ids, function(x) paste(x, '_abundance', sep='')))
colnames(txi_gene$length) <- unlist(lapply(sample_ids, function(x) paste(x, '_length', sep='')))

if (save==TRUE) {
    log_print(paste(Sys.time(), 'Saving txi files...'))
    if (!dir.exists(file.path(wd, 'data', 'malaria', 'txi'))) {
        dir.create(file.path(wd, 'data', 'malaria', 'txi'), recursive=TRUE)
    }
    write.table(txi_gene$counts,
                file.path(wd, 'data', 'malaria', 'txi', 'txi_counts.csv'),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')
    write.table(txi_gene$abundance,
                file.path(wd, 'data', 'malaria', 'txi', 'txi_abundances.csv'),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')
    # write.table(txi_gene$length,
    #             file.path(wd, 'data', 'malaria', 'txi', 'txi_lengths.csv'),
    #             quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')
}

myDGEList <- DGEList(txi_gene$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)


# ----------------------------------------------------------------------
# EDA
# from malarialab_solution.R


# Plot unfiltered, un-normalized log2 cpm
log_print(paste(Sys.time(), 'Plotting unfiltered, un-normalized log2 cpm...'))

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sample_ids)
log2.cpm.df.pivot <- pivot_longer(
    log2.cpm.df, # dataframe to be pivoted
    cols = -1, # column names to be stored as a SINGLE variable
    names_to = "samples", # name of that new variable (column)
    values_to = "expression"
) # name of new variable (column) storing all the values (data)

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
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'malaria', 'eda', 'cpm_unfiltered_nonnormalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# Plot filtered, non-normalized counts per million
log_print(paste(Sys.time(), 'Plotting filtered, non-normalized cpm...'))

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>1)>=2 #user defined
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sample_ids)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = -1, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
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

if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'malaria', 'eda', 'cpm_filtered_nonnormalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# Plot filtered, normalized cpm
log_print(paste(Sys.time(), 'Plotting filtered, normalized cpm...'))

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sample_ids)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = -1, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

if (save==TRUE) {
    log_print(paste(Sys.time(), 'Saving filtered, normalized cpm...'))
    write.table(
        log2.cpm.filtered.norm,
        file.path(wd, 'data', 'malaria', 'filtered_normalized_cpm.csv'),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep=','
    )
}

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
    ggsave(file.path(wd, 'figures', 'malaria', 'eda', 'tpm_filtered_normalized.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# Plot all violins together
log_print(paste(Sys.time(), 'Plotting cowplot...'))
cowplot::plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)
if (save==TRUE) {
    ggsave(file.path(wd, 'figures', 'malaria', 'eda', 'cowplot.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()


# ----------------------------------------------------------------------
# Code graveyard

# txi_gene = join_many_csv(
#     paths,
#     index_cols=c('target_id'),
#     value_cols=c('length', 'eff_length', 'est_counts', 'tpm'),
#     sep='\t'
# )

# # Pivot to long format
# tpm_cols = filter_list_for_match(colnames(txi_gene), pattern=c('tpm'))
# tpms <- melt(
#     txi_gene[, c('target_id', tpm_cols)],
# )
# colnames(tpms) = c('target_id', 'sample_id', 'tpm')
