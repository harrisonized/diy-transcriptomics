wd = dirname(this.path::here())  # wd = '~/github/diy-transcriptomics'
source(file.path(wd, 'R', 'utils.R'))
# library(tidyverse)  # too broad
library('readr')  # read_delim
library('dplyr')
library('stringr')  # could not find function "str_detect"
library('ggplot2')
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
log <- log_open(paste0("lemis_eda-", strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))



# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading data...'))
lemis <- read_delim(file.path(wd, "data", "lemis", "lemis_cleaned.tsv"))
# lemis <- read.csv(
#     file.path(wd, "data", "lemis", "lemis_cleaned.tsv"),
#     header=TRUE, sep='\t', check.names=FALSE
# )  # slower

# inspect
# glimpse(lemis)
# colnames(lemis)
# typeof(lemis)  # list


# ----------------------------------------------------------------------
# EDA


# Identify the most common (by ‘quantity’) live mammal taken from the wild for import into the US.
# Building on your analysis above, produce a plot showing live mammals (use ‘generic_name’)
# imported for the purposes of science/research. (tip: use geom_col() in ggplot for this).
# Feel free to play around with different themes to make your plot more exciting.

log_print(paste(Sys.time(), 'Generating Figure 1...'))
agg_tbl <- lemis %>%
    filter(source=='Specimens taken from the wild',
           taxa=='mammal',
           description=='Live specimens (live animals or plants)',
           purpose=='Scientific'  # additional filter
           ) %>%
    group_by(genus, species, specific_name, generic_name) %>%
    summarize(total_quantity = sum(quantity)) %>%
    arrange(desc(total_quantity))
# A tibble: 45 × 5
# Groups:   genus, species, specific_name [45]
#    genus          species      specific_name  generic_name total_quantity
#    <chr>          <chr>        <chr>          <chr>                 <int>
#  1 Rousettus      aegyptiacus  EGYPTIAN       ROUSETTE                510
#  2 Macaca         fascicularis CRAB-EATING    MACAQUE                 478
#  3 Chlorocebus    aethiops     GRIVET         MONKEY                  182
#  4 Heterocephalus glaber       NAKED          MOLE-RAT                155
#  5 Aotus          nancymaae    NIGHT          MONKEY                  140
#  6 Martes         pennanti     NA             FISHER                   96
#  7 Peromyscus     maniculatus  NORTH AMERICAN DEERMOUSE                74
#  8 Bison          bison        AMERICAN       BISON                    73
#  9 Scotinomys     sp.          SINGING        MOUSE                    72
# 10 Cephalophus    maxwelli     MAXWELL'S      DUIKER                   65
fig <- ggplot(
    agg_tbl[1:10, ], aes(x = reorder(generic_name, -total_quantity), y = total_quantity)) +
    geom_bar(stat = "identity") + 
    xlab("generic name") +
    ylab("total_quantity") +
    ggtitle("Live mammals Imported for Science") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'lemis', 'figure_1-imported_live_mammals.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# Identify the countries from which we import the most macaques (again, a simple plot will suffice).
log_print(paste(Sys.time(), 'Generating Figure 2...'))
agg_tbl <- lemis %>%
    filter(generic_name=='MACAQUE') %>%
    group_by(country_origin) %>%
    summarize(total_quantity = sum(quantity)) %>%
    arrange(desc(total_quantity))
fig <- ggplot(agg_tbl[1:10, ], aes(x = reorder(country_origin, -total_quantity), y = total_quantity)) +
    xlab("country_origin") +
    geom_bar(stat = "identity") + 
    ggtitle("Countries we import the most macaques") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'lemis', 'figure_2-countries_most_macaques.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# Using the same approach as above, create a plot showing the countries from which we import live bats.
log_print(paste(Sys.time(), 'Generating Figure 3...'))
agg_tbl <- lemis %>%
    filter(generic_name=='BAT',
           description=='Live specimens (live animals or plants)') %>%
    group_by(country_origin) %>%
    summarize(total_quantity = sum(quantity)) %>%
    arrange(desc(total_quantity))
fig <- ggplot(agg_tbl[1:10, ], aes(x = reorder(country_origin, -total_quantity), y = total_quantity)) +
    xlab("country_origin") +
    geom_bar(stat = "identity") + 
    ggtitle("Countries we import the most live bats") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'lemis', 'figure_3-countries_most_bats.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# For what purposes do we import bats?
log_print(paste(Sys.time(), 'Generating Figure 4...'))
agg_tbl <- lemis %>%
    filter(generic_name=='BAT',
           description=='Live specimens (live animals or plants)') %>%
    group_by(purpose) %>%
    summarize(total_quantity = sum(quantity)) %>%
    arrange(desc(total_quantity))
fig <- ggplot(agg_tbl[1:10, ], aes(x = reorder(purpose, -total_quantity), y = total_quantity)) +
    xlab("purpose") +
    geom_bar(stat = "identity") + 
    ggtitle("Purpose that we import the most live bats") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'lemis', 'figure_4-countries_most_bats_purpose.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# How does the type of bat (use ‘specific_name’) imported differ between countries
# (hint: use facet_wrap in your ggplot code)?
log_print(paste(Sys.time(), 'Generating Figure 5...'))
agg_tbl <- lemis %>%
    filter(generic_name=='BAT',
           description=='Live specimens (live animals or plants)') %>%
    group_by(country_origin, specific_name) %>%
    summarize(total_quantity = sum(quantity)) %>%
    arrange(desc(total_quantity))
fig <- ggplot(agg_tbl, aes(fill=specific_name, y=total_quantity, x=reorder(country_origin, -total_quantity))) + 
    xlab("purpose") +
    geom_bar(position="stack", stat="identity") + 
    ggtitle("Types of imported bats from different countries") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'lemis', 'figure_5-imported_bats.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# You are alerted to a concerning new viral disease of humans that is believed to originate from Fruit bats
# (though the exact type of fruit bat is not clear). Identify the US cities that would would be most likely
# to be exposed to such a virus from the import of live fruit bats.
log_print(paste(Sys.time(), 'Generating Figure 6...'))
agg_tbl <- lemis %>%
    filter(generic_name=='BAT',
           description=='Live specimens (live animals or plants)',
           str_detect(.$specific_name, "FRUIT")
        ) %>%
    group_by(port, country_origin, specific_name) %>%
    summarize(total_quantity = sum(quantity)) %>%
    arrange(desc(total_quantity))

fig <- ggplot(agg_tbl, aes(fill=specific_name, y=total_quantity, x=reorder(port, -total_quantity))) + 
    xlab("purpose") +
    geom_bar(position="stack", stat="identity") + 
    ggtitle("Types of imported fruit bats from different countries") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'lemis', 'figure_6-imported_fruit_bats.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# A recent case of Anthrax in NYC was traced back to a contaminated Wildebeest hide that was stretched and
# used to make a traditional drum. Through which port(s) did this animal product most likely enter the country?
# (note: this actually happens).
log_print(paste(Sys.time(), 'Generating Figure 7...'))
agg_tbl <- lemis %>%
    filter(generic_name=='WILDEBEEST',
           taxa=='mammal',
           str_detect(.$description, "skin")
        ) %>%
    group_by(port, country_origin, specific_name) %>%
    summarize(total_quantity = sum(quantity)) %>%
    arrange(desc(total_quantity))
fig <- ggplot(agg_tbl, aes(fill=specific_name, y=total_quantity, x=reorder(port, -total_quantity))) + 
    xlab("port") +
    geom_bar(position="stack", stat="identity") + 
    ggtitle("Types of imported wildebeests at different ports") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
if (!troubleshooting) {
    ggsave(file.path(wd, 'figures', 'lemis', 'figure_7-wildebeest-skin.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()


# ----------------------------------------------------------------------
# Code graveyard

# Identify the most common (by ‘quantity’) live mammal taken from the wild for import into the US.
# agg_tbl <- lemis %>%
#     filter(source=='Specimens taken from the wild',
#            taxa=='mammal',
#            description=='Live specimens (live animals or plants)',
#            ) %>%
#     group_by(genus, species, specific_name, generic_name) %>%
#     summarize(total_quantity = sum(quantity)) %>%
#     arrange(desc(total_quantity))

# A tibble: 390 × 5
# Groups:   genus, species, specific_name [389]
#    genus      species      specific_name  generic_name total_quantity
#    <chr>      <chr>        <chr>          <chr>                 <int>
#  1 Ursus      americanus   AMERICAN BLACK BEAR                  46087
#  2 Lepus      americanus   SNOWSHOE       RABBIT                17090
#  3 Macaca     fascicularis CRAB-EATING    MACAQUE               14409
#  4 Odocoileus virginianus  WHITE-TAILED   DEER                   4332
#  5 Dasyprocta sp.          NA             AGOUTI                 3133
#  6 NA         NA           NON-CITES      MAMMALS                3018
#  7 Grammomys  sp.          NA             RAT                    3000
#  8 Cynictis   sp.          NA             MONGOOSE               2852
#  9 Cervus     elaphus      NA             ELK                    2639
# 10 Petaurus   breviceps    SUGAR          GLIDER                 2460
# … with 380 more rows
# ℹ Use `print(n = ...)` to see more rows


# Identify the most expensive (by ‘value’) shipment of live mammals to enter the US.
# agg_tbl <- lemis %>%
#     filter(taxa=='mammal',
#            description=='Live specimens (live animals or plants)',
#            ) %>%
#     arrange(desc(value))
# agg_tbl[which.max(agg_tbl$value),]
# Giant Panda

# How does the answer above compare with the most expensive shipment of any kind (live or not)?
# lemis[which.max(lemis$value),]
# Asian Elephant Ivory pieces (not manufactured, includes scraps)
