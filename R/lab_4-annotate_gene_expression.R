library('EnsDb.Hsapiens.v86')
library('BSgenome.Mfuro.UCSC.musFur1')
library('biomaRt')

# set working directory
wd = '/home/harrisonized/github/diy-transcriptomics'


# ----------------------------------------------------------------------
# Utils

list_files <- function(dir_path, ext=NULL, recursive = TRUE) {
    all_files = list.files(dir_path, recursive = recursive, full.name=TRUE)

    if (!is.null(ext)) {
        # See: https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
        return (all_files[tools::file_ext(all_files)==ext])
    } else {
        return (all_files)
    }
}

filter_list_for_match <- function(items, pattern) {
    # filter
    for (i in 1:length(pattern)){
        items <- lapply(items, grep, pattern=pattern[[i]], value=TRUE)
    }
    return (unlist(items[!sapply(items, identical, character(0))]))  # remove character(0)
}


# ----------------------------------------------------------------------
# Load transcripts

Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- transcripts(BSgenome.Mfuro.UCSC.musFur1, columns=c("tx_id", "gene_name"))  # doesn't work, because of the BSgenome


# ----------------------------------------------------------------------
# Get from biomart

# myMart <- useMart(biomart="ENSEMBL_MART_ENSEMBL")
# available.datasets <- listDatasets(myMart)

# get ferret genome
ferret.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "mpfuro_gene_ensembl")
ferret.attributes <- listAttributes(ferret.anno)
Tx.ferret <- getBM(
    attributes=c('ensembl_transcript_id_version', 'external_gene_name'),
    mart = ferret.anno
)  # query the biomaRt database

Tx.ferret <- tidyr::as_tibble(Tx.ferret)  # convert to tibble
Tx.ferret <- dplyr::rename(
    Tx.ferret, target_id = ensembl_transcript_id_version, 
    gene_name = external_gene_name
)  # rename columns


# annotate
files = filter_list_for_match(list_files(file.path(wd, "data", "lab_2", "mapped_reads")), 'h5')
Txi_gene <- tximport::tximport(
    files,  # load data from mapped_reads
    type = "kallisto", 
    tx2gene = Tx,
    txOut = TRUE,  # How does the result change if this =FALSE vs =TRUE?  # if this is false, it doesn't work
    countsFromAbundance = "lengthScaledTPM",
    ignoreTxVersion = TRUE
)


# ----------------------------------------------------------------------
# Get ferret promoter genes

# Help
# ?getSequence

# This is just doing a SQL query
# IFIT2, OAS2, IRF1, IFNAR1, and MX1
sequences = getSequence(
    id = c("IFIT2", "OAS2", "IRF1", "IFNAR1", "MX1"),
    type = "external_gene_name",  # or "hgnc_symbol" and "uniprot_gn_symbol" are missing MX1
    seqType = "coding_gene_flank",
    upstream = 1000,
    mart = ferret.anno
)

# Error in .processResults(postRes, mart = mart, hostURLsep = sep, fullXmlQuery = fullXmlQuery,  : 
# Query ERROR: caught BioMart::Exception::Database: 
# Error during query execution: You have an error in your SQL syntax; 
# check the manual that corresponds to your MySQL server version for the right syntax to use near
# 'AND main.transcript_id_1064_key=mpfuro_gene_ensembl__exon_transcript__dm.transcr' at line 1
