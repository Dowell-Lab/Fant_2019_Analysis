.libPaths( c( .libPaths(), "/Users/zama8258/R") )
library(tidyverse)

## Set up parallel processing
library(parallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)

dir <- "/scratch/Users/zama8258/DMSO2_3_FIX.BedGraph"
print(dir)
out <- str_c(dir, ".flat", sep="")
print(out)

## Initialize Dataset
dat <- read_tsv(file=dir, col_names = c("chrom", "start", "end", "reads"))

dat$coord <- mcmapply(seq,dat$start,dat$end,SIMPLIFY=FALSE)
flatdat <- dat %>%
    unnest(coord) %>%
            select(-start,-end)

write_tsv(flatdat, out, append=FALSE)
