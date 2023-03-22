# ask is to go through the Step1 through Step6 scripts and generate some figures

library(ggplot2)
library(tximport)
library(reshape)

wd = '/home/harrisonized/github/diy-transcriptomics'


# Functions
# # join_many_csv
# # filter_list_for_match


join_many_csv <- function(filepaths, index_cols, value_cols, ext='csv', recursive=TRUE, sep=',') {
    filenames = c(tools::file_path_sans_ext(basename(dirname(paths))))  # get the foldername

    
    # read dfs and left join on index_cols
    df_list <- lapply(filepaths, read.csv, sep=sep)

    # Warning: column names ‘count.x’, ‘count.y’ are duplicated in the result
    # See: https://stackoverflow.com/questions/38603668/suppress-any-emission-of-a-particular-warning-message
    withCallingHandlers({
        all_reads <- Reduce(
            function(...) merge(..., by=index_cols),
            lapply(df_list, "[", c(index_cols, value_cols))
        )
    }, warning = function(w) {
        # print(conditionMessage(w))
        if (startsWith(conditionMessage(w), "column names")) {
            invokeRestart( "muffleWarning" )
        }
    })
    
    # rename columns
    colnames(all_reads) = c(
        index_cols,  # index_cols
        as.list(outer(value_cols, filenames, paste, sep='-'))  # suffix value_cols with filename
    )
    return(all_reads)
}


filter_list_for_match <- function(items, pattern) {
    # filter
    for (i in 1:length(pattern)){
        items <- lapply(items, grep, pattern=pattern[[i]], value=TRUE)
    }
    return (unlist(items[!sapply(items, identical, character(0))]))  # remove character(0)
}



# ----------------------------------------------------------------------
# Read Data

paths = file.path(
    list.dirs(
        path = file.path(wd, 'data', 'lab_9', 'malaria'),
        full.names = TRUE,
        recursive = FALSE
    ),
    "abundance.tsv"
)

all(file.exists(path))

# wtf man
# txi_gene <- tximport(
#     paths,
#     type = "kallisto",
#     txOut = TRUE, #How does the result change if this =FALSE vs =TRUE?
#     countsFromAbundance = "lengthScaledTPM",
#     ignoreTxVersion = TRUE
# )


txi_gene = join_many_csv(
    paths,
    index_cols=c('target_id'),
    value_cols=c('length', 'eff_length', 'est_counts', 'tpm'),
    sep='\t'
)
tpm_cols = filter_list_for_match(colnames(txi_gene), pattern=c('tpm'))


# need to pivot first...
tpms <- melt(
    txi_gene[, c('target_id', tpm_cols)],
    id = c("sample_id")
)
colnames(tpms) = c('target_id', 'sample_id', 'tpm')



# ----------------------------------------------------------------------
# Plot the unnormalized, unfiltered TPM first

p1 <- ggplot(tpms) +
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
p1

ggsave(
    "/home/harrisonized/github/diy-transcriptomics/figures/lab_9/tpm.png",
    scale = 1,
    width = 800,
    height = 600,
    units = "px",
    dpi = 100,
    plot = p1
)

ggsave(
    "/home/harrisonized/github/diy-transcriptomics/figures/lab_9/tpm-1.5.png",
    scale = 1.5,
    width = 800,
    height = 600,
    units = "px",
    dpi = 150,
    plot = p1
)


# ----------------------------------------------------------------------
# Need to make volcano plot
# Need to do PCA
# Need to plot heatmap

