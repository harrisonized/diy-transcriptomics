## Source: Step1_TxImport.R
## This script imports the kalliso outputs from the leishmania dataset
## using 'EnsDb.Hsapiens.v86' to annotate each row of the resulting dataframe.
## No outputs, run this in RStudio

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
suppressMessages(library('EnsDb.Hsapiens.v86'))
# suppressMessages(library('BSgenome.Mfuro.UCSC.musFur1'))
suppressMessages(library('GenomicFeatures'))
library('tibble')
library('optparse')
library('logr')
source(file.path(wd, 'R', 'functions', 'utils.R'))  # list_files, filter_list_for_match


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default="data/leishmania/mapped_reads",
                metavar="data/leishmania/mapped_reads", type="character",
                help="set the input directory"),

    make_option(c("-f", "--filename"), default="abundance.h5",
                metavar="abundance.h5", type="character",
                help="Choose: 'abundance.h5' or 'abundance.tsv'. '.h5' is faster for importing data"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt['troubleshooting'][[1]]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("leishmania_eda-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data


# define data source
files = filter_list_for_match(
    list_files(file.path(wd, opt['input-dir'][[1]])),
    opt['filename'][[1]]  # file_ext
)
# Example file: "data/schistosoma/mapped_reads/F12h_LE_1/abundance.tsv"
#   target_id    length eff_length est_counts   tpm
#   <chr>         <dbl>      <dbl>      <dbl> <dbl>
# 1 Smp_186980.1    402      239.       171   58.5 
# 2 Smp_197050.1    186       50.3       59   95.8 
# 3 Smp_160500.1   5316     5152.       244    3.87


# get human annotations
tx2gene_obj <- as_tibble(
    GenomicFeatures::transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
)  
tx2gene_obj <- dplyr::rename(tx2gene_obj, target_id = tx_id)  # change first column name to 'target_id'
tx2gene_obj <- dplyr::select(tx2gene_obj, "target_id", "gene_name")  # transcript ID needs to be the first column in the dataframe


# Read abundance files into a txi object
# see: https://github.com/mikelove/tximport/blob/devel/R/tximport.R
# names(txi)
# [1] "abundance"  "counts"  "infReps"  "length"  "countsFromAbundance
# txi$abundance = tpm
# txi$counts != counts, they did something
# txi$infReps: inferential replicates, not sure what this is
# txi$length = eff_length
# txi$countsFromAbundance = "lengthScaledTPM"
sample_txi <- tximport::tximport(
    files,  # load data from mapped_reads
    type = "kallisto", 
    tx2gene = tx2gene_obj,  # use EnsDb.Hsapiens.v86 to annotate the data
    txOut = TRUE,  # if this is false, it doesn't work
    countsFromAbundance = "lengthScaledTPM",
    ignoreTxVersion = TRUE
)


# ----------------------------------------------------------------------
# Peform EDA here



end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
