## Source: Step1_TxImport.R
## This script imports the kalliso outputs from the leishmania dataset
## using 'EnsDb.Hsapiens.v86' to annotate each row of the resulting dataframe.
## No outputs, run this in RStudio

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
suppressMessages(library('EnsDb.Hsapiens.v86'))
# suppressMessages(library('BSgenome.Mfuro.UCSC.musFur1'))
suppressMessages(library('GenomicFeatures'))
library('optparse')
library('logr')
source(file.path(wd, 'R', 'utils.R'))  # list_files, filter_list_for_match


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default="data/leishmania/mapped_reads",
                metavar="data/leishmania/mapped_reads", type="character",
                help="set the input directory"),

    make_option(c("-f", "--file-ext"), default="h5",
                metavar="h5", type="character",
                help="Choose: 'h5' or 'tsv'. 'h5' is faster for importing data"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="disable if for troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt['troubleshooting'][[1]]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("leishmania_annotate_gene_expression-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Import Kallisto data

# define data sources
files = filter_list_for_match(
    list_files(file.path(wd, opt['input-dir'][[1]])),
    opt['file-ext'][[1]]  # file_ext
)
tx2gene_obj <- GenomicFeatures::transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))

# annotate, returns a list of dataframes
sample_txi <- tximport::tximport(
    files,  # load data from mapped_reads
    type = "kallisto", 
    tx2gene = tx2gene_obj,
    txOut = TRUE,  # if this is false, it doesn't work
    countsFromAbundance = "lengthScaledTPM",
    ignoreTxVersion = TRUE
)

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
