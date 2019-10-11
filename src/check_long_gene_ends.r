## Script to run DESeq2
## Maintainer: Zachary Maas <zama8258@colorado.edu>

#######################
## Preliminary Setup ##
#######################

library("tidyverse")
library("ggthemes")
library("DESeq2")

## Change this for different experiments...
setwd('/home/zach/dowell_lab/pausing_meta_analysis/out/normalization')
num_reps <- 2

## To filter out only the genes we've used and convert to common ID
idConvert <- read_delim("refseq_to_common_id.txt",
                        col_names = c("rowname", "common"), delim = " ")
## Experimental condition
condition <- factor(c(rep("treated", num_reps), rep("untreated", num_reps)), levels=c("untreated", "treated"))

## We have to fix the counts table by removing the first row, hence "counts_fix.txt"
full_seqdata <- read.delim("counts_fix.txt", stringsAsFactors = FALSE, header = TRUE,
                           row.names = 1)
## Then, we filter to only include the resolved isoforms
full_df <- as_tibble(rownames_to_column(full_seqdata))

## For DROSOPHILA
if (!is.null(full_df$rowname)) {
    full_df$common <- full_df$rowname
    full_df$rowname <- NULL
}

full_dt <- data.frame(full_df)
rownames(full_dt) <- full_dt$common
full_dt$common <- NULL
full_dt <- as_tibble(full_dt)

full_dt <- full_dt %>% mutate(c_mean = (rowMeans(full_dt[c('C413_1_S3_R1_001.sorted.bam',
                                                           'C413_2_S4_R1_001.sorted.bam')])),
                              p_mean = (rowMeans(full_dt[c('PO_1_S1_R1_001.sorted.bam',
                                                           'PO_2_S2_R1_001.sorted.bam')])))

## Characterize the underlying normalized read distribution
ggplot(data = full_dt) +
    geom_density(aes(x = c_mean, color="Control")) +
    geom_density(aes(x = p_mean, color="Treatment")) +
    labs(title = "Distribution of Coverage in the 120kb -> -500polyA Region",
         x = "Coverage",
         y = "Proportion",
         color = "Sample") +
    theme_tufte() +
    scale_x_log10() +
    annotation_logticks() +
    ggsave("dist.pdf", width = 10, height = 5, device = "pdf")

## Use linear regression to check for baseline skew between samples
reg <- function(df) {
    lin_mod <- lm(c_mean ~ p_mean, data=df)
    slope <- summary(lin_mod)$coefficients[2]
    rsqr <- summary(lin_mod)$r.squared
    return(slope)
}
slope <- reg(full_dt)
msg <- paste0("Treatment/Control Slope: ", signif(slope, 5))

## NOTE: This step is not required, it's just for checking the quality
## of results.

## Perform monte carlo subsampling to verify those results. Zhang, P.
## (1993). Model Selection Via Muiltfold Cross Validation. Ann. says
## that using N^2 will give us close to optimal results. That will
## take too long on my laptop, so we'll just do 100k trials, which is
## about 10% of that value.
sample_slopes <- replicate(100000, reg(sample_frac(full_dt, 0.25)))
write(sample_slopes, "slopes.txt")
max_distance <- max(sample_slopes)
min_distance <- min(sample_slopes)
write(paste0("Subsampled 100000x at 25% Original Size\n",
             "Max: ", max_distance, "\n",
             "Min: ", min_distance), "subsample_info.txt")

## Generate a distribution plot for the cross validation
ggplot() +
    geom_density(aes(x = sample_slopes)) +
    labs(title = "Distribution of Regression Slopes from Cross Validation",
         x = "Slope",
         y = "Proportion") +
    theme_tufte() +
    xlim(0, 2) +
    ggsave("crossval_dist.pdf", width = 10, height = 5, device = "pdf")

## Generate a plot showing the linear regression.
ggplot(data = full_dt) +
    geom_point(aes(x = c_mean, y = p_mean)) +
    geom_smooth(method = 'lm', aes(x = c_mean, y = p_mean)) +
    geom_abline(slope=1, color = "red") +
    annotate("text", x = 10, y = 1000, label = msg) +
    labs(title = "Regression For Normalization Based on +120kb to -0.5kb on Genes > 120kb",
         x = "Control Counts",
         y = "Treatment Counts") +
    theme_tufte() +
    scale_x_log10() + scale_y_log10() +
    ggsave("corr.pdf", width = 10, height = 5, device = "pdf")
