## Generate figures...

load("/scratch/Users/zama8258/pi_test.Rda")

dd <- as.tibble(pi) %>% filter(GeneID != "")

dd$GeneID = as.integer(as.character(dd$GeneID))

tsd  <- as.tibble(ts) 
colnames(tsd)[colnames(tsd) == "tx_id"] <- "GeneID"

prd2k <- read_delim("/scratch/Users/zama8258/genefilter_test/pause_ratios.data",
                    col_names=c('tx_name', 'strand', 'Pause_2k'),
                    delim=" ")
prd5k <- read_delim("/scratch/Users/zama8258/genefilter_test/pause_ratios_5k.data",
                    col_names=c('tx_name', 'strand', 'Pause_5k'),
                    delim=" ")

ddt <- as.tibble(
    merge(x=dd, y=tsd[, c("tx_name", "GeneID")],
          by.x="GeneID", by.y="GeneID", all.x = TRUE))
ddt <- ddt[, c("GeneID", "tx_name", "Fisher", "Pause")]
ddt$tx_name <- substr(ddt$tx_name,1,nchar(ddt$tx_name)-2)
ddt <- ddt %>% arrange(tx_name)
prd <- prd %>% arrange(tx_name)
ddt <- left_join(x=ddt, y=prd2k[, c("tx_name", "Pause_2k")], by="tx_name")
ddt <- left_join(x=ddt, y=prd5k[, c("tx_name", "Pause_5k")], by="tx_name")
ddt <- na.omit(ddt)

ddt  <- ddt %>% filter(Fisher < 0.05)

ddf <- ddt %>% filter(Pause_2k < 200, Pause_5k < 200)

library(ggthemes)
library(ggsci)
library(wesanderson)

ggplot(data = ddt, mapping=aes(x=Pause_2k, y=Pause_5k), alpha=1/10) +
    geom_jitter() + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() +
    labs(x = "2kb Window",
         y= "5kb Window",
         title= "Comparison of Pausing Index Methods") 
ggsave("/scratch/Users/zama8258/pause_comparison.png", plot = last_plot(), device = "png")

ggplot(data = ddf, mapping=aes(x=Pause_2k, y=Pause_5k), alpha=1/10) +
    geom_point() + geom_abline(aes(intercept=0,slope=1)) +
    theme_tufte() + xlim(0,200) + ylim(0,200) +
    labs(x = "2kb Window",
         y= "5kb Window",
         title= "Comparison of Pausing Index Methods (PI < 200)") 
ggsave("/scratch/Users/zama8258/pause_comparison_filt.png", plot = last_plot(), device = "png")
