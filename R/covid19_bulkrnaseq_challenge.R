## See: https://diytranscriptomics.com/lab/lab-10

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
# library(tidyverse)  # too broad
library('readr')
library('textshape')
# library('ggplot2')
library('optparse')
library('logr')


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-h", "--human-data"),
                default="data/covid19_bulkrnaseq_challenge/GSE147507_RawReadCounts_Human.tsv",
                metavar="data/covid19_bulkrnaseq_challenge/GSE147507_RawReadCounts_Human.tsv",
                type="character",
                help="path/to/human_data.tsv"),

    make_option(c("-a", "--animal-data"),
                default="data/covid19_bulkrnaseq_challenge/GSE147507_RawReadCounts_Ferret.tsv",
                metavar="data/covid19_bulkrnaseq_challenge/GSE147507_RawReadCounts_Ferret.tsv",
                type="character",
                help="path/to/animal_data.tsv"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt['troubleshooting'][[1]]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("covid19_bulkrnaseq_challenge-", strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

human_data <- read_tsv(file.path(wd, opt['human-data'][[1]]))
human_data <- as.matrix(column_to_rownames(human_data, "...1"))

animal_data <- read_tsv(file.path(wd, opt['ferret-data'][[1]]))
animal_data <- as.matrix(column_to_rownames(animal_data, "...1"))


# ----------------------------------------------------------------------
# Peform EDA here


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
