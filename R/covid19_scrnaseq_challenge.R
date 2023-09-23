## See: https://diytranscriptomics.com/lab/lab-13

wd = dirname(this.path::here())  # wd = '~/github/R/diy-transcriptomics'
suppressMessages(library('Seurat'))
library('Matrix')
library('tibble')
suppressMessages(library('dplyr'))
library('rjson')
library('R2HTML')
library('readr')
suppressMessages(library('textshape'))
suppressMessages(library('scales'))
library('ggplot2')
suppressMessages(library('DT'))
library('optparse')
library('logr')
source(file.path(wd, 'R', 'functions', 'utils.R'))  # load_rdata


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"),
                default="data/covid19_scrnaseq_challenge/immuneHealth_integrated.gz",
                metavar="data/covid19_scrnaseq_challenge/immuneHealth_integrated.gz",
                type="character", help="path to seurat object"),

    make_option(c("-m", "--metadata"),
                default="data/covid19_scrnaseq_challenge/ImmuneHealth_COVID_studyDesign.txt",
                metavar="data/covid19_scrnaseq_challenge/ImmuneHealth_COVID_studyDesign.txt",
                type="character", help="path/to/studyDesign.txt"),

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
log <- log_open(paste0("covid19_scrnaseq_challenge-", strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

seurat_obj <- load_rdata(file.path(wd, opt['input-file'][[1]]))
study_design <- read_tsv(file.path(wd, opt['metadata'][[1]]))


# ----------------------------------------------------------------------
# Peform EDA here


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
