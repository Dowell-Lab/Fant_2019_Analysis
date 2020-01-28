library('tidyverse')

sites <- c(31827806, 31827822, 31827878, 31827969)
process <- function(filename) {
    dat <- read_delim(filename,
                      col_names=c('chr', 'start', 'end', 'count'),
               delim='\t') %>%
        mutate(selected = ifelse(start %in% sites, "SITE", "NOTSITE")) %>%
        group_by(selected)
    return(dat)
}

ct <- process('reads_ct.bed')
ko <- process('reads_ko.bed')

ct_sum <- ct %>% summarise(sum(count))
ko_sum <- ko %>% summarise(sum(count))

ct_vals <- ct_sum$`sum(count)`
ko_vals <- ko_sum$`sum(count)`

ct_vals[2] / ct_vals[1]
ko_vals[2] / ko_vals[1]
