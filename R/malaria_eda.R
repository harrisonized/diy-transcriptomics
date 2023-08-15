# Lab 9
# Ask is to go through the Step_1 through Step_6 scripts and generate some figures


wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
source(file.path(wd, 'R', 'utils.R'))
# library(tidyverse) # too broad
library('ggplot2')
library('tximport')
library('reshape')
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
# Read Data

paths = filter_list_for_match(
    list_files(file.path(wd, 'data', 'malaria')),
    'abundance.tsv'
)
txi_gene = join_many_csv(
    paths,
    index_cols=c('target_id'),
    value_cols=c('length', 'eff_length', 'est_counts', 'tpm'),
    sep='\t'
)


# ----------------------------------------------------------------------
# Plot the unnormalized, unfiltered TPM first

# Pivot to long format
tpm_cols = filter_list_for_match(colnames(txi_gene), pattern=c('tpm'))
tpms <- melt(
    txi_gene[, c('target_id', tpm_cols)],
)
colnames(tpms) = c('target_id', 'sample_id', 'tpm')

fig <- ggplot(tpms) +
    aes(x=sample_id, y=tpm, fill=sample_id) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    stat_summary(
        fun = "median", 
        geom = "point", 
        shape = 95, 
        size = 10, 
        color = "black", 
        show.legend = FALSE) +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) +
    ylim(0, 250) +
    labs(y="tpm", x = "sample_id",
       title="Transcripts Per Million",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time()))
if (save==TRUE) {
    ggsave(
        file.path(wd, 'figures', 'malaria', 'tpm.png'),
        height = 600, width = 800, scale = 1.5, units = "px", dpi = 150,
        plot = fig
    )
}


# ----------------------------------------------------------------------
# Need to make volcano plot
# Need to do PCA
# Need to plot heatmap

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
