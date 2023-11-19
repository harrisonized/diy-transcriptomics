## Adapted from: DIY_scRNAseq.R
## Filters empty drops and creates some QC plots

wd = dirname(this.path::here())  # wd = '~/github/R/diy-transcriptomics'
suppressMessages(library('Seurat'))
library('Matrix')
suppressMessages(library('DropletUtils'))
library('tibble')
suppressMessages(library('dplyr'))
library('rjson')
library('R2HTML')
library('readr')
suppressMessages(library('textshape'))
library('scales')
library('ggplot2')
suppressMessages(library('DT'))
library('optparse')
library('logr')
source(file.path(wd, 'R', 'functions', 'scrnaseq_qc_plots.R'))  # bc_rank_plot, print_HTML
source(file.path(wd, 'R', 'functions', 'utils.R'))  # read_10x


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"),
                default='data/covid19_scrnaseq/pbmc_1k_v3_scRNAseq_processed',
                metavar='data/covid19_scrnaseq/pbmc_1k_v3_scRNAseq_processed',
                type="character",
                help="cellranger directory containing barcodes.tsv, genes.tsv, and matrix.mtx"),
 
    make_option(c("-o", "--output-dir"), default="figures/covid19_scrnaseq",
                metavar="figures/covid19_scrnaseq", type="character",
                help="set the output directory for the figures"),

    make_option(c("-f", "--fdr-threshold"), default=0.05, metavar="0.05",
                type="double", help="set the max fdr threshold"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

# create required directories
if (!troubleshooting) {
    if (!dir.exists(file.path(wd, opt[['output-dir']]))) {
        dir.create(file.path(wd, opt[['output-dir']]))
    }
    if (!dir.exists(file.path(wd, opt[['input-dir']], 'qc_report'))) {
        dir.create(file.path(wd, opt[['input-dir']], 'qc_report'))
    }
}

# Start Log
start_time = Sys.time()
log <- log_open(paste0("covid19_scrnaseq_qc",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading data...'))

raw_mtx <- read_10x(
    file.path(wd, opt[['input-dir']], 'counts_unfiltered', 'cellranger')
)

# load run info from JSON files produced by Kb
kb_stats <- c(
    fromJSON(file = file.path(wd, opt[['input-dir']], 'inspect.json')),
    fromJSON(file = file.path(wd, opt[['input-dir']], 'run_info.json'))
)


# ----------------------------------------------------------------------
# Filter empty drops

log_print(paste(Sys.time(), 'Filtering empty drops...'))

# use DropletUtils package to get probability that each barcode is a cell
# this takes a minute
drop_stats <- DropletUtils::emptyDrops(raw_mtx)

# set threshold probability for calling a cell
filt_mtx <- raw_mtx[,(
    (drop_stats[['FDR']] <= opt[['fdr-threshold']]) &
    (is.na(drop_stats[['FDR']])==FALSE)
)]
filt_barcodes = data.frame('V1'=colnames(filt_mtx))  # columns
filt_genes = data.frame('X2'=rownames(filt_mtx))  # columns

# write drop_stats filtered results
if (!troubleshooting) {
    write10xCounts(
        file.path(wd, opt[['input-dir']], 'counts_filtered', 'cellranger'),
        filt_mtx,
        gene.symbol = filt_genes[['X2']],
        overwrite=TRUE  # warning: this will empty your directory!
    )
}


# ----------------------------------------------------------------------
# Generate QA report

# this report will contain some useful metrics as well as the traditional
# log-transformed UMI rank plot (a.k.a. 'waterfall' plot)
# plot was first described in the Drop-seq paper: - Macosko et al. 2015, DOI:10.1016/j.cell.2015.05.002
# this plot has two important points that we will try to identify:
# 1. 'knee point' - is the point where the signed curvature is minimized. 
# This corresponds to a transition between a distinct subset of barcodes
# with large totals and the majority of barcodes with smaller totals
# 2. 'inflection point' - is the point on the curve where the first derivative is minimized. 
# This corresponds to the point past which cells cannot reliably be distinguished from background

log_print(paste(Sys.time(), 'Generating QC report...'))

# determine chemistry version
tech <- grep('10X(.*)', unlist(strsplit(kb_stats[['call']], '\\s')), value=TRUE) 

# make a nice/simple table that summarizes that stats
seq_stats <- data.frame(
    # get sequencing/alignment stats 
    stat = c('Sequencing technology', 'Number of reads processed',
             '% reads pseudoaligned', '% reads on whitelist'),
    value = prettyNum(
        c(tech, kb_stats[['n_processed']], kb_stats[['p_pseudoaligned']], 
          round(kb_stats[['percentageReadsOnWhitelist']],2)),
    big.mark = ',')
)

if (!troubleshooting) {
    write.table(
        seq_stats,
        file.path(wd, opt[['input-dir']], 'qc_report', 'seq_stats.csv'),
        row.names = FALSE,
        sep = ','
    )
}

# calculate cell stats and save to df
pct_counts_in_cells <- round((sum(filt_mtx)/sum(raw_mtx))*100, 2) 
med_counts_cell <- median(colSums(filt_mtx))
med_genes_cell <- median(apply(filt_mtx, 2, function(x) sum(x >= 1)))
tot_genes_detected <- sum(rowSums(filt_mtx)>=1)
cell_stats <- data.frame(
    stat = c('Estimated number of cells', '% counts in cells', 
             'Median counts per cell', 'Median genes per cell',
             'Total genes detected'), 
    value = prettyNum(
        c(ncol(filt_mtx), pct_counts_in_cells, med_counts_cell,
          med_genes_cell, tot_genes_detected),
    big.mark = ',')
)

if (!troubleshooting) {
    write.table(
        cell_stats,
        file.path(wd, opt[['input-dir']], 'qc_report', 'cell_stats.csv'),
        row.names = FALSE,
        sep = ','
    )
}

# create barcode rank plot png
if (!troubleshooting) {
    bc_rank_plot(
        stats = barcodeRanks(raw_mtx),
        raw_cells = raw_barcodes, filt_cells = filt_barcodes,
        save = file.path(wd, opt[['output-dir']], 'barcode_rank.png')
    )
}

# output a HTML summary of the run
# currently broken because references in the function are off
# print_HTML(
#     seq_stats = seq_stats, cell_stats = cell_stats,
#     dir = file.path(wd, 'data', 'covid19_scrnaseq', 'pbmc_1k_v3_scRNAseq_processed',
#               'counts_unfiltered'),
#     sample_id = NULL
# )

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
