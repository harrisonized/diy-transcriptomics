## Common use functions

# Functions
# # list_files
# # filter_list_for_match
# # join_many_csv


# list all files in a directory
list_files <- function(dir_path, ext=NULL, recursive = TRUE) {
    all_files = list.files(dir_path, recursive = recursive, full.name=TRUE)

    if (!is.null(ext)) {
        # See: https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
        return (all_files[tools::file_ext(all_files)==ext])
    } else {
        return (all_files)
    }
}


# python equivalent: [string for string in items if re.match(pattern, string)]
filter_list_for_match <- function(items, pattern) {
    for (i in 1:length(pattern)){
        items <- lapply(items, grep, pattern=pattern[[i]], value=TRUE)
    }
    return (unlist(items[!sapply(items, identical, character(0))]))  # remove character(0)
}


# loads an RData file, and returns it
# Without this function, base R uses the filename as the variable name
# See: https://stackoverflow.com/questions/5577221/can-i-load-a-saved-r-object-into-a-new-object-name
load_rdata <- function(filepath){
    load(filepath)
    return( get(basename(filepath)) )
}


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
