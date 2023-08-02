## Common use functions

# Functions
# # list_files
# # filter_list_for_match


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
