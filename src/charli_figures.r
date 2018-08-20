suppressMessages(library("tidyverse"))
library(ggthemes)
## suppressMessages(library(GenomicFeatures))
## suppressMessages(library(GenomicRanges))
## suppressMessages(library(rtracklayer))

## Generate figures...

c_1_201 <- read_delim("C413_1_S3_R1_001.trim.rpkm.bedGraph_pause_ratios_201.data",
                      col_names=c('tx_name', 'strand', 'c_1_201'), delim=" ") %>%
    subset(select=-strand)
c_2_201 <- read_delim("C413_2_S4_R1_001.trim.rpkm.bedGraph_pause_ratios_201.data",
                      col_names=c('tx_name', 'strand', 'c_2_201'), delim=" ") %>%
    subset(select=-strand)
p_1_201 <- read_delim("PO_1_S1_R1_001.trim.rpkm.bedGraph_pause_ratios_201.data",
                      col_names=c('tx_name', 'strand', 'p_1_201'), delim=" ") %>%
    subset(select=-strand)
p_2_201 <- read_delim("PO_2_S2_R1_001.trim.rpkm.bedGraph_pause_ratios_201.data",
                      col_names=c('tx_name', 'strand', 'p_2_201'), delim=" ") %>%
    subset(select=-strand)

c_1_301 <- read_delim("C413_1_S3_R1_001.trim.rpkm.bedGraph_pause_ratios_301.data",
                      col_names=c('tx_name', 'strand', 'c_1_301'), delim=" ") %>%
    subset(select=-strand)
c_2_301 <- read_delim("C413_2_S4_R1_001.trim.rpkm.bedGraph_pause_ratios_301.data",
                      col_names=c('tx_name', 'strand', 'c_2_301'), delim=" ") %>%
    subset(select=-strand)
p_1_301 <- read_delim("PO_1_S1_R1_001.trim.rpkm.bedGraph_pause_ratios_301.data",
                      col_names=c('tx_name', 'strand', 'p_1_301'), delim=" ") %>%
    subset(select=-strand)
p_2_301 <- read_delim("PO_2_S2_R1_001.trim.rpkm.bedGraph_pause_ratios_301.data",
                      col_names=c('tx_name', 'strand', 'p_2_301'), delim=" ") %>%
    subset(select=-strand)

c_1_5000 <- read_delim("C413_1_S3_R1_001.trim.rpkm.bedGraph_pause_ratios_5000.data",
                       col_names=c('tx_name', 'strand', 'c_1_5000'), delim=" ") %>%
    subset(select=-strand)
c_2_5000 <- read_delim("C413_2_S4_R1_001.trim.rpkm.bedGraph_pause_ratios_5000.data",
                       col_names=c('tx_name', 'strand', 'c_2_5000'), delim=" ") %>%
    subset(select=-strand)
p_1_5000 <- read_delim("PO_1_S1_R1_001.trim.rpkm.bedGraph_pause_ratios_5000.data",
                       col_names=c('tx_name', 'strand', 'p_1_5000'), delim=" ") %>%
    subset(select=-strand)
p_2_5000 <- read_delim("PO_2_S2_R1_001.trim.rpkm.bedGraph_pause_ratios_5000.data",
                       col_names=c('tx_name', 'strand', 'p_2_5000'), delim=" ") %>%
    subset(select=-strand)

ddt_201 <- list(c_1_201, c_2_201, p_1_201, p_2_201) %>% reduce(left_join, by="tx_name")
ddt_201$c_mean <- rowMeans(ddt_201[c('c_1_201', 'c_2_201')], na.rm=TRUE)
ddt_201$p_mean <- rowMeans(ddt_201[c('p_1_201', 'p_2_201')], na.rm=TRUE)
ddt_201 <- ddt_201 %>% subset(select=c(tx_name, c_mean, p_mean))

ddt_301 <- list(c_1_301, c_2_301, p_1_301, p_2_301) %>% reduce(left_join, by="tx_name")
ddt_301$c_mean <- rowMeans(ddt_301[c('c_1_301', 'c_2_301')], na.rm=TRUE)
ddt_301$p_mean <- rowMeans(ddt_301[c('p_1_301', 'p_2_301')], na.rm=TRUE)
ddt_301 <- ddt_301 %>% subset(select=c(tx_name, c_mean, p_mean))

ddt_5000 <- list(c_1_5000, c_2_5000, p_1_5000, p_2_5000) %>% reduce(left_join, by="tx_name")
ddt_5000$c_mean <- rowMeans(ddt_5000[c('c_1_5000', 'c_2_5000')], na.rm=TRUE)
ddt_5000$p_mean <- rowMeans(ddt_5000[c('p_1_5000', 'p_2_5000')], na.rm=TRUE)
ddt_5000 <- ddt_5000 %>% subset(select=c(tx_name, c_mean, p_mean))

ddt_201 <- ddt_201 %>% filter(c_mean > 0.01, p_mean > 0.01)

## 201 scatter
ggplot(data = ddt_201, mapping=aes(x=p_mean, y=c_mean), alpha=1/10) +
    geom_point() + geom_abline(aes(intercept=0,slope=1)) +
    ## theme_tufte() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/201_scatter.png"),
    plot = last_plot(), device = "png")

## 201 ecdf
ggplot(data = ddt_201) +
    stat_ecdf(aes(p_mean, color='red')) + stat_ecdf(aes(c_mean, color='blue')) +
    ## theme_tufte() +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/201_ecdf.png"),
    plot = last_plot(), device = "png")

## 301 scatter
ggplot(data = ddt_301, mapping=aes(x=p_mean, y=c_mean), alpha=1/10) +
    geom_point() + geom_abline(aes(intercept=0,slope=1)) +
    ## theme_tufte() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_scatter.png"),
    plot = last_plot(), device = "png")

## 301 ecdf
ggplot(data = ddt_301, alpha=1/10) +
    stat_ecdf(aes(p_mean)) + stat_ecdf(aes(c_mean)) +
    ## theme_tufte() +
    labs(x = "Pausing Index",
         y= "Cumulative Distribution",
         title= "Cumulative Distribution") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/301_ecdf.png"),
    plot = last_plot(), device = "png")

## 5000 scatter
ggplot(data = ddt_5000, mapping=aes(x=p_mean, y=c_mean), alpha=1/10) +
    geom_bin2d(bins=100) + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    labs(x = "Control",
         y= "Treatment",
         title= "Change in Pausing Index") +
    scale_x_log10() + scale_y_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/5000_scatter.png"),
    plot = last_plot(), device = "png")

## 5000 ecdf
ggplot(data = ddt_5000, alpha=1/10) +
    stat_ecdf(aes(p_mean, color='Control')) +
    stat_ecdf(aes(c_mean, color='Treatment')) +
    theme_tufte() +
    theme(legend.title=element_blank()) +
    labs(x = "Pausing Index",
         y= "Cumulative Density",
         title= "Cumulative Density Distribution of Pausing Indices") +
    scale_x_log10()
ggsave(
    str_c("/scratch/Users/zama8258/pause_output/5000_ecdf.png"),
    plot = last_plot(), device = "png")
