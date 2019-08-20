### explore_lengths.r --- examine lengths for analysis
##
## Filename: explore_lengths.r
## Created: Mon Jun  3 15:45:12 2019 (-0600)
##
######################################################################
##
### Commentary:
##
## An exploration of gene lengths.
##
######################################################################
##
### Code:

library('tidyverse')
library('ggthemes')

## setwd('/scratch/Users/zama8258/processed_nascent/metagene/')
setwd('/home/zach/dowell_lab/pausing_meta_analysis/out/length_experiment/')

geneimport <- function(file) {
    q <- read_delim(file,
                    col_names = c('chr', 'start', 'end', 'name', 'x', 'strand'),
               delim = '\t') %>% mutate(len = abs(end - start))
    return(q)
}

q1 <- geneimport("quartile1Genes.bed")
q2 <- geneimport("quartile2Genes.bed")
q3 <- geneimport("quartile3Genes.bed")
q4 <- geneimport("quartile4Genes.bed")

summary(q1$len)
summary(q2$len)
summary(q3$len)
summary(q4$len)

ggplot() +
    theme_tufte() +
    geom_density(dat = q1, aes(x = len, color = 'Quartile 1')) +
    geom_density(dat = q2, aes(x = len, color = 'Quartile 2')) +
    geom_density(dat = q3, aes(x = len, color = 'Quartile 3')) +
    geom_density(dat = q4, aes(x = len, color = 'Quartile 4')) +
    labs(title = 'Distribution of Gene Lengths',
         x = 'Gene Length (bp)', y = 'Frequency',
         color = 'Quartile') +
    xlim(c(0, 200000)) +
    ggsave('quartile_distribution_plot.pdf',
           width = 10, height = 5)

## Part 2 of Analysis

r <- read_delim('sorted.fpkm',
                col_names = c('chr', 'start', 'end', 'name', 'x', 'strand', 'expression', 'geneid'),
                delim = '\t') %>% mutate(len = abs(end - start))

ggplot() +
    theme_tufte() +
    geom_point(dat = r, aes(x = len, y = expression)) +
    labs(title = 'Gene Length vs Normalized Expression',
         x = 'Gene Length (bp)', y = 'Expression',
         color = 'Quartile') +
    scale_x_log10() +
    scale_y_log10() +
    ggsave('length_vs_expression_plot.pdf',
           width = 10, height = 5)

######################################################################
### explore_lengths.r ends here
