## Source: lab4_solution.R
## Example for using getBM() and getSequences() to query the bioMart database
## This is a standalone script

wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
library('tibble')
suppressMessages(library('biomaRt'))
suppressMessages(library('EnsDb.Hsapiens.v86'))
library('optparse')
library('logr')


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-m", "--biomart-dataset"), default="mpfuro_gene_ensembl",
                metavar="mpfuro_gene_ensembl", type="character",
                help="the class used 'mpfuro_gene_ensembl'"),

    make_option(c("-o", "--output-dir"), default="data/biomart",
                metavar="data/biomart", type="character",
                help="set the output directory"),

    make_option(c("-g", "--genes"), default=c("IFIT2", "OAS2", "IRF1", "IFNAR1", "MX1"),
                metavar='c("IFIT2", "OAS2", "IRF1", "IFNAR1", "MX1")', type="list",
                help="genes of interest"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="disable if for troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt['troubleshooting'][[1]]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("query_biomart-", strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Get reference transcripts from biomart

log_print(paste(Sys.time(), 'Getting reference dataset...'))

# query the biomaRt database
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = opt['biomart-dataset'][[1]])
ref_tx <- getBM(
    attributes=c('ensembl_transcript_id_version', 'external_gene_name'),
    mart = mart
)
ref_tx <- tidyr::as_tibble(ref_tx)  # convert to tibble
ref_tx <- dplyr::rename(ref_tx, 
                        target_id = ensembl_transcript_id_version,
                        gene_name = external_gene_name)  # rename columns

# save
if (troubleshooting==FALSE){
    write.table(
        ref_tx,
        file = file.path(wd, opt['output-dir'][[1]],  # dirpath
                         paste0(opt['biomart-dataset'][[1]], '.csv')),  # filename
        row.names = FALSE,
        sep = ','
    )
}

# examine the mart object
# attributes <- listAttributes(mart)  


# ----------------------------------------------------------------------
# Query promoter genes

log_print(paste(Sys.time(), 'Querying promoter sequences...'))

sequences = getSequence(
    id = unlist(opt['genes'][[1]]),
    type = "external_gene_name",  # or "hgnc_symbol" and "uniprot_gn_symbol" are missing MX1
    seqType = "coding_gene_flank",
    upstream = 1000,
    mart = mart
)
sequences <- as_tibble(rev(sequences))

# save
if (troubleshooting==FALSE){
    write.table(
        sequences,
        file = file.path(wd, opt['output-dir'][[1]],  # dirpath
                         paste0('promoter_sequences', '.csv')),  # filename
        row.names = FALSE,
        sep = ','
    )
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()