suppressMessages(library("tidyverse"))
library(ggthemes)
library(ggsci)

## Generate figures...
## You only need to adjust the filenames here...
c_1_301 <- read_delim("C413_1_S3_R1_001.trim.rpkm.bedGraph_pause_ratios_301.data",
                     col_names=c('tx_name', 'strand', 'c_1_301', 'coverage_c1_301'), delim=" ") %>%
    subset(select=-strand)
c_2_301 <- read_delim("C413_2_S4_R1_001.trim.rpkm.bedGraph_pause_ratios_301.data",
                      col_names=c('tx_name', 'strand', 'c_2_301', 'coverage_c2_301'), delim=" ") %>%
    subset(select=-strand)
p_1_301 <- read_delim("PO_1_S1_R1_001.trim.rpkm.bedGraph_pause_ratios_301.data",
                      col_names=c('tx_name', 'strand', 'p_1_301', 'coverage_p1_301'), delim=" ") %>%
    subset(select=-strand)
p_2_301 <- read_delim("PO_2_S2_R1_001.trim.rpkm.bedGraph_pause_ratios_301.data",
                      col_names=c('tx_name', 'strand', 'p_2_301', 'coverage_p2_301'), delim=" ") %>%
    subset(select=-strand)

ddt_301 <- list(c_1_301, c_2_301, p_1_301, p_2_301) %>% reduce(left_join, by="tx_name")
ddt_301$c_pause_mean <- rowMeans(ddt_301[c('c_1_301', 'c_2_301')], na.rm=TRUE)
ddt_301$p_pause_mean <- rowMeans(ddt_301[c('p_1_301', 'p_2_301')], na.rm=TRUE)
ddt_301$c_coverage_mean <- rowMeans(ddt_301[c('coverage_c1_301', 'coverage_c2_301')], na.rm=FALSE)
ddt_301$p_coverage_mean <- rowMeans(ddt_301[c('coverage_p1_301', 'coverage_p2_301')], na.rm=FALSE)

ddt_301 <- subset(ddt_301, !duplicated(tx_name))
ddt_301 <- ddt_301 %>% mutate(pause_diff = log2(c_pause_mean) - log2(p_pause_mean))
ddt_301$coverage_mean <- rowMeans(ddt_301[c('c_coverage_mean', 'p_coverage_mean')], na.rm=TRUE)

## 301 scatter
ggplot(data = ddt_301, mapping=aes(x=p_pause_mean, y=c_pause_mean), alpha=1/10) +
geom_bin2d(bins=75) + geom_abline(aes(intercept=0,slope=1)) +
theme_tufte() +
scale_fill_gsea() +
labs(x = "Control",
     y= "Treatment",
     title= "Change in Pausing Index") +
scale_x_log10() + scale_y_log10()
ggsave(
    str_c("301_scatter.png"),
    plot = last_plot(), device = "png")

## 301 ecdf
ggplot(data = ddt_301) +
stat_ecdf(aes(p_pause_mean, color='Control')) + stat_ecdf(aes(c_pause_mean, color='Treatment')) +
theme_tufte() +
scale_color_npg() +
theme(legend.title=element_blank()) +
labs(x = "Pausing Index",
     y= "Cumulative Distribution",
     title= "Cumulative Distribution") +
scale_x_log10()
ggsave(
    str_c("301_ecdf.png"),
    plot = last_plot(), device = "png")

fddt_301 <- filter(ddt_301, pause_diff > -6e-5)

## MA-ish plot
ggplot(data = ddt_301, mapping=aes(x=coverage_mean, y=pause_diff), alpha=1/10) +
geom_bin2d(bins=75) + geom_hline(yintercept=0) +
scale_fill_gsea() +
theme_tufte() +
labs(x = "Normalized Coverage Level",
     y= "Change in Pause Index",
     title= "MA (log2(treatment) - log2(control)) vs Coverage Level (log10)") +
scale_x_log10() +
## scale_y_log10()
ggsave(
    str_c("301_MA.png"),
    plot = last_plot(), device = "png")
