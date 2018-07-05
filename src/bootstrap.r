## Bootstrapping script to install all necessary packages

.libPaths( c( .libPaths(), "/Users/zama8258/R/") )
.libPaths()
message("Established libpaths")

## Install Bioconductor Packages
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite("GenomicFeatures")
biocLite("rtracklayer")
biocLite("groHMM")
message("Installed Bioconductor ")

install.packages(c("tidyverse", "argparse", "ggthemes"))

library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(groHMM)
