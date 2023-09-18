# See: https://diytranscriptomics.com/lab/lab-10

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
# library(tidyverse) # too broad
library('readr')
library('textshape')
# library('ggplot2')
library('optparse')
library('logr')


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt['troubleshooting'][[1]]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("covid19_eda-", strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

human_covid_data <- read_tsv(file.path(wd, "data", "covid19_challenge", "GSE147507_RawReadCounts_Human.tsv"))
human_covid_data <- as.matrix(column_to_rownames(human_covid_data, "...1"))

ferret_covid_data <- read_tsv(file.path(wd, "data", "covid19_challenge", "GSE147507_RawReadCounts_Human.tsv"))
ferret_covid_data <- as.matrix(column_to_rownames(ferret_covid_data, "...1"))

# ----------------------------------------------------------------------
# Peform EDA here



end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
