## Preliminary Pausing Analysis
## Maintainter: Zachary Maas <zama8258@colorado.edu>
## Licensed under GPLv3

################################################################################
########################## Install and Load Libraries ##########################
################################################################################

.libPaths( c( .libPaths(), "/Users/zama8258/R/") )
.libPaths()
message("Established libpaths")

## TODO - Extract to installation sctipt
## source("https://bioconductor.org/biocLite.R")
## biocLite("Rsubread")

## install.packages("wesanderson")

## General utilities to make R behave better
library(tidyverse)

## Gives us access to FeatureCount
library("Rsubread")

## Fast rolling means
library(RcppRoll)

## Set up parallel processing
library(parallel)

## Generate a cluster for parallel processing
no_cores <- detectCores()
## no_cores <- 16
print(str_c("Processing with ", no_cores, " cores"))

## cl <- makeCluster(no_cores)

## Assign our bedfile from the parsed arguments
## TODO - Implement argparse
analysis_bedfile <- "/scratch/Shares/public/nascentdb/processedv2.0/bedgraphs/SRR2084576.tri.BedGraph"

################################################################################
####################### Extract Regions Proximal to TSS ########################
################################################################################

message("Conducting featureCount Analysis")
## Featurecount analysis
## TODO - Need to pre-convert from a bed using tempfile?
fc <- featureCounts(
    files="/scratch/Shares/public/nascentdb/processedv2.0/bams/SRR2084577.trimmed.bam",
    annot.inbuilt="hg19",
    GTF.featureType="exon",
    strandSpecific=0,
    ## Should increase this...
    nthreads=no_cores)

annot <- as_data_frame(fc$annotation)

## Take only the 5' most site...
annot$Chr <- gsub(";.*","",annot$Chr)
annot$Start <- as.numeric(gsub(";.*","",annot$Start))
annot$End <- as.numeric(gsub(";.*","",annot$End))
annot$Strand <- gsub(";.*","",annot$Strand)

## TODO Make this user adjustable, not just 2500
annot$StartExt <- annot$Start-2500
annot$EndExt <- annot$End+2500

## Extract our regions of interest
interest_regions <- select(annot, Chr, StartExt, EndExt)

## Filter out non-numerical chromosomes
## TODO add back in x/y
interest_regions_chrs <- interest_regions %>%
    filter(grepl('chr\\d+', Chr)) %>% 
    mutate_if(is.numeric, as.integer)

################################################################################
############## Sort Our Regions of Interest for Bedtools Intersect #############
################################################################################

message("Filtering Regions of Interest using Bedtools")

## Tempfile for our regions of interest
tmp_interest <- tempfile()
write_tsv(interest_regions_chrs, tmp_interest,
          append=FALSE, col_names=FALSE) 

## Tempfile for sorting our bedfile
tmp_sorted <- tempfile()
sort_cmd <- str_c("module load bedtools/2.25.0 &&", "bedtools sort",
                 "-i", tmp_interest,
                 ">", tmp_sorted,
                 sep=" ")

## Execute command
try(system(sort_cmd))

################################################################################
################ Intersect Regions of Interest with User Bedfile ###############
################################################################################

## Tempfile for when we run bedtools intersect
tmp_intersected <- tempfile()

## Command to run our bedtools intersect on our region of interest
inter_cmd <- str_c("module load bedtools/2.25.0 &&", "bedtools intersect",
                  "-sorted", "-wa",
                  "-a", analysis_bedfile,
                  "-b", tmp_interest,
                  ">", tmp_intersected,
                  sep=" ")

try(system(inter_cmd))

################################################################################
######################## Parse Newly Intersected Bedfile #######################
################################################################################

dat <- read_tsv(tmp_intersected, col_names = c("chrom", "start", "end", "reads"))
## dat <- read_tsv(analysis_bedfile, col_names = c("chrom", "start", "end", "reads"))

## Split by strand
pos_dat <- filter(dat, dat$reads > 0)
neg_dat <- filter(dat, dat$reads < 0)

################################################################################
############################### Helper Functions ###############################
################################################################################

genomeFlatten <- function(bedData) {
    rng_pre <- (mcmapply(`:`, bedData$start, bedData$end))

    ## We need this to make R behave well when we have a small number
    ## of rows (1/2), since it tries to coerce those rows in the wrong
    ## direction...
    ## TODO - Document Better
    if (!is.null(ncol(rng_pre))) {
        rng_pre <- split(rng_pre, rep(1:ncol(rng_pre), each = nrow(rng_pre)))
        }     
    bedData["rng"] <- list(rng_pre)

    if (nrow(bedData) == 0) {
        message("[MSG] Cannot flatten - no data.")
        return(NA)
    }
    
    ## TODO - Can we speed this up at all?
    flatdat <- bedData %>% bind_rows(bedData) %>%    # make larger sample data
        mutate_if(is.list, simplify_all) %>%    # flatten each list element internally 
        unnest() %>%    # expand
        select(-start, -end)

    flatdat <- flatdat %>%
        complete(rng = full_seq(flatdat$rng, 1), chrom,
                 fill=list(reads=0)) %>%
        distinct()

    return(flatdat)
}

genomeFlatten2 <- compiler::cmpfun(genomeFlatten)

################################################################################
########################## Functions to calculate TSS ##########################
################################################################################

findSenseTSS <- function(data, refseq_tss_coord, range, win_size, step){

    coord_start <- refseq_tss_coord - range
    coord_end <- refseq_tss_coord + range

    ## Filter our data to only be within ±`range` of the tss coordinate
    data_in_tss_window <- filter(data,
                                data$start > coord_start - 15000,
                                data$end < coord_end + 15000)

    ## Flatten our data for the rolling mean
    flatdat <- genomeFlatten2(data_in_tss_window)
    
    ## Return early if we have problems finding our TSS. This will
    ## happen if it turns out that we have no reads within the window
    ## that we are looking for the TSS in.
    if (is.na(flatdat) || nrow(flatdat) == 0) {
        ## TODO Check NA vs 0
        print(str_c("[MSG] Could not find TSS for: ", refseq_tss_coord))
        return(0)
    }

    ## If we don't have enough reads to do the full rolling mean, just
    ## take the 5' most read value
    if (nrow(flatdat) < 50) {
        return(min(flatdat$rng))
    }
    
    ## Find Highest Mean, convert to TSS coordinate
    filt_rollmean <- roll_mean(flatdat$reads, by=step, n=win_size)

    if (length(filt_rollmean) == 0) {
        ## TODO How else can we handle this?
        print(str_c("[MSG] Insufficent data to find TSS for: ", refseq_tss_coord))
        return(0)
    }
    
    ## Find the number of windows from 0 that get us to the maximum mean,
    ## treated as the TSS.  NOTE - this returns a 1 indexed-value
    max_win_pos <- grep(max(filt_rollmean), filt_rollmean)


    ## Find the tss coordinate
    tss_coord_is <- min(data_in_tss_window$start) + (step * (max_win_pos-1))

    return(first(tss_coord_is))
}

findAntiSenseTSS <- function(data, refseq_tss_coord, range, win_size, step){

    coord_start <- refseq_tss_coord - range
    coord_end <- refseq_tss_coord + range

    ## Filter our data to only be within ±`range` of the tss coordinate
    data_in_tss_window <- filter(data,
                                data$start > coord_start - 15000,
                                data$end < coord_end + 15000)

    ## Flatten our data for the rolling mean
    flatdat <- genomeFlatten2(data_in_tss_window)
    
    ## Return early if we have problems finding our TSS. This will
    ## happen if it turns out that we have no reads within the window
    ## that we are looking for the TSS in.
    if (is.na(flatdat) || nrow(flatdat) == 0) {
        ## TODO Check NA vs 0
        print(str_c("[MSG] Could not find TSS for: ", refseq_tss_coord))
        return(0)
    }

    ## If we don't have enough reads to do the full rolling mean, just
    ## take the 5' most read value
    if (nrow(flatdat) < 50) {
        return(min(flatdat$rng))
    }
    
    ## Find Highest Mean, convert to TSS coordinate
    filt_rollmean <- roll_mean(flatdat$reads, by=step, n=win_size)

    if (length(filt_rollmean) == 0) {
        ## TODO How else can we handle this?
        print(str_c("[MSG] Insufficent data to find TSS for: ", refseq_tss_coord))
        return(0)
    }
    
    ## Find the number of windows from 0 that get us to the maximum mean,
    ## treated as the TSS.  NOTE - this returns a 1 indexed-value
    max_win_pos <- grep(max(filt_rollmean), filt_rollmean)


    ## Find the tss coordinate
    tss_coord_is <- min(data_in_tss_window$start) + (step * (max_win_pos-1))

    return(first(tss_coord_is))
}

################################################################################
####################### Functions to Determine Gene Body #######################
################################################################################

findSenseGeneBody <- function(data, start_coord,
                             pause_window_size,
                             threshold, cutoff) {

    ## Establish variables
    scan_start <- start_coord + pause_window_size
    ## TODO - Find a better constant...
    coord_start <- start_coord - 15000
    coord_end <- start_coord + 15000

    ## Preliminary Filtering to make sure we keep data
    data_in_tss_window <- filter(data,
                                data$start >= coord_start,
                                data$end < coord_end)    
    

    ## If we don't have any reads, then our gene body doesn't exist
    if (nrow(data_in_tss_window) == 0) {
        return(NA)
    }

    ## Flatten our data by individual coordinates
    flatdat <- genomeFlatten2(data_in_tss_window)

    ## Secondary filtering to only the range we care about
    data_in_tss_window <- filter(data,
                                data$start >= scan_start,
                                data$end < scan_start + cutoff)    


    ## TODO - Fix wrangling on this
    if (nrow(flatdat) == 0) {
        ## 
        print(str_c("[MSG] Could not determine gene body for: ", start_coord))
        return(NA)
    }

    ## Iterate, checking for if we make the cutoff.
    ## TODO - Mary Allen's Comment About Length...
    currCoord <- start_coord
    totDist <- 0
    tmpDist <- 0
    currCount <- 0
    cutoff_cond <- FALSE
    while (cutoff_cond == FALSE) {
        ## Find the count of our item and add it to the count variable
        tmp_pos <- which(currCoord == flatdat$rng)
        ## tmp_pos <- Position(p(`==`, currCoord), flatdat$rng)
        if (length(tmp_pos) == 0) {
            tmp_sum <- 0
        } else {
            tmp_sum <- flatdat$reads[tmp_pos]
        }
        
        currCount <- currCount + first(tmp_sum)

        ## Increment total and current lookup variables
        tmpDist <- tmpDist + 1
        currCoord <- currCoord + 1

        ## If we don't have enough reads and go outside the cutoff,
        ## we're done
        if (currCount < threshold && tmpDist >= cutoff) {
            cutoff_cond <- TRUE
        }
        ## If we do have enough reads, add the distance to totDist and
        ## reset tmpDist
        if (currCount >= threshold) {
            totDist <- totDist + tmpDist
            tmpDist <- 0
            currCount <- 0
        }
    }

    ## Return how far we were able to go from the TSS
return(totDist)
}

findAntiSenseGeneBody <- function(data, start_coord,
                                 pause_window_size,
                             threshold, cutoff) {

    ## Establish variables
    scan_start <- start_coord - pause_window_size
    ## TODO - Find a better constant...
    coord_start <- start_coord - 15000
    coord_end <- start_coord + 15000

    ## Preliminary Filtering to make sure we keep data
    ## print(data)
    ## data_in_tss_window <- filter(data,
    ##                             data$start >= coord_end,
    ##                             data$end < coord_start)    

    ## If we don't have any reads, then our gene body doesn't exist
    if (nrow(data) == 0) {
        return(NA)
    }

    ## Flatten our data by individual coordinates
    flatdat <- genomeFlatten2(data)
    print(flatdat)

    ## Secondary filtering to only the range we care about
    data_in_tss_window <- filter(data,
                                data$start <= scan_start,
                                data$end > scan_start - cutoff)    


    ## TODO - Fix wrangling on this
    if (nrow(flatdat) == 0) {
        ## 
        print(str_c("[MSG] Could not determine gene body for: ", start_coord))
        return(NA)
    }

    ## Iterate, checking for if we make the cutoff.
    ## TODO - Mary Allen's Comment About Length...
    currCoord <- start_coord
    totDist <- 0
    tmpDist <- 0
    currCount <- 0
    cutoff_cond <- FALSE
    while (cutoff_cond == FALSE) {
        ## Find the count of our item and add it to the count variable
        tmp_pos <- which(currCoord == flatdat$rng)
        ## tmp_pos <- Position(p(`==`, currCoord), flatdat$rng)
        if (length(tmp_pos) == 0) {
            tmp_sum <- 0
        } else {
            tmp_sum <- flatdat$reads[tmp_pos]
        }
        
        currCount <- currCount + first(tmp_sum)

        ## Increment total and current lookup variables
        tmpDist <- tmpDist - 1
        currCoord <- currCoord - 1

        ## If we don't have enough reads and go outside the cutoff,
        ## we're done
        if (currCount > threshold && tmpDist >= cutoff) {
            cutoff_cond <- TRUE
        }
        ## If we do have enough reads, add the distance to totDist and
        ## reset tmpDist
        if (currCount >= threshold) {
            totDist <- totDist - tmpDist
            tmpDist <- 0
            currCount <- 0
        }
    }

    ## Return how far we were able to go from the TSS
    return(totDist)
}

################################################################################
######################## Parsing the Reference Sequence ########################
################################################################################

## Reference Sequence TSS Parsing
## ## TODO - Parameterize Filename
refseq <- read_tsv("/scratch/Users/zama8258/NCBI_RefSeq_UCSC_RefSeq_hg19.bed",
                  col_names = c("chrom", "start", "end", "name", "score",
                                "strand", "thickstart", "thickend", "itemRGB",
                                "blockCount", "blockSizes", "blockStarts"))

## ## Separate by strand and take only what we need.
pos_tss <- refseq %>% filter(strand == '+')
pos_tss <- pos_tss %>% subset(select=c("chrom", "start","end"))
neg_tss <- refseq %>% filter(strand == '-')
neg_tss <- neg_tss %>% subset(select=c("chrom", "start","end"))

tss_test <- pos_tss[!duplicated(pos_tss$end),]
tss_test <- tss_test %>%
    filter(grepl('chr\\d+|chrY|chrX', chrom))

################################################################################
##################### Function to Calculate Pausing Ratio ######################
################################################################################

calcSensePauseRatio <- function(data, refseq_tss_coord, tss_range,
                               tss_win_size, tss_calc_step,
                          body_threshold, body_cutoff) {

    message(str_c("[COORD] ", refseq_tss_coord))
    ## Flatten data for later processing...
    ## TODO - Change 15000 to user parameter
    data_we_care_about <- filter(data,
                                data$start >= refseq_tss_coord - tss_win_size - 15000,
                                data$end < refseq_tss_coord + 15000)
    flatdat <- genomeFlatten2(data_we_care_about)

    ## We need to suppress warnings here, since this comparison will
    ## throw a warning for non-empty dataframes
    suppressWarnings(
        if (is.na(flatdat)) {
            message(str_c("[MSG] No reads found proximal to: ", refseq_tss_coord))
            return(c(refseq_tss_coord, 0, 0, NA))
        }
    )
    ## Find our TSS start coordinate 
    start_coord <- findSenseTSS(data_we_care_about, refseq_tss_coord,
                               tss_range, tss_win_size, tss_calc_step)

    ## Bail out early if we fail
    if (start_coord == 0) {
        return(c(refseq_tss_coord, 0, 0, NA))
    }
    
    ## Filter our data so we can find a sum for TSS reads
    data_in_tss <- filter(flatdat,
                         flatdat$rng >= start_coord,
                         flatdat$rng < start_coord + tss_win_size)
    tss_reads <- sum(data_in_tss$reads)

    ## Find Gene Body Length
    body_length <- findSenseGeneBody(data_we_care_about, start_coord,
                                    tss_win_size, body_threshold, body_cutoff) 

    ## Filter again to find a sum for gene body reads
    data_in_gene_body <- filter(flatdat,
                               flatdat$rng >= start_coord + tss_win_size,
                               flatdat$rng < start_coord +
                               tss_win_size + body_length)
    gene_body_reads <- sum(data_in_gene_body$reads)

    ## Calculate and yield pausing index
    if (tss_reads != 0 && gene_body_reads != 0) {
        pausing_index <- tss_reads / gene_body_reads
    } else {
        pausing_index <- NA
    }

    return(
        list(refseq_tss_coord, tss_reads, gene_body_reads, pausing_index))
}

calcSensePauseRatio2 <- compiler::cmpfun(calcSensePauseRatio)

calcAntiSensePauseRatio <- function(data, refseq_tss_coord, tss_range,
                                   tss_win_size, tss_calc_step,
                               body_threshold, body_cutoff) {

    message(str_c("[COORD] ", refseq_tss_coord))
    ## Flatten data for later processing...
    ## TODO - Change 15000 to user parameter
    data_we_care_about <- filter(data,
                                data$start <= refseq_tss_coord - tss_win_size + 15000,
                                data$end > refseq_tss_coord - 15000)
    ## print(data_we_care_about)
    flatdat <- genomeFlatten2(data_we_care_about)

    ## We need to suppress warnings here, since this comparison will
    ## throw a warning for non-empty dataframes
    suppressWarnings(
        if (is.na(flatdat)) {
            message(str_c("[MSG] No reads found proximal to: ", refseq_tss_coord))
            return(c(refseq_tss_coord, 0, 0, NA))
        }
    )
    ## Find our TSS start coordinate 
    start_coord <- findAntiSenseTSS(data_we_care_about, refseq_tss_coord,
                                   tss_range, tss_win_size, tss_calc_step)

    ## Bail out early if we fail
    if (start_coord == 0) {
        return(c(refseq_tss_coord, 0, 0, NA))
    }
    
    ## Filter our data so we can find a sum for TSS reads
    data_in_tss <- filter(flatdat,
                         flatdat$rng <= start_coord,
                         flatdat$rng > start_coord - tss_win_size)
    tss_reads <- sum(data_in_tss$reads)

    ## Find Gene Body Length
    body_length <- findAntiSenseGeneBody(data_we_care_about, start_coord,
                                        tss_win_size, body_threshold, body_cutoff) 

    print(body_length)
    ## Filter again to find a sum for gene body reads
    data_in_gene_body <- filter(flatdat,
                               flatdat$rng >= start_coord + tss_win_size,
                               flatdat$rng < start_coord +
                               tss_win_size + body_length)
    gene_body_reads <- sum(data_in_gene_body$reads)

    ## Calculate and yield pausing index
    if (tss_reads != 0 && gene_body_reads != 0) {
        pausing_index <- tss_reads / gene_body_reads
    } else {
        pausing_index <- NA
    }

    ## print(list(refseq_tss_coord, tss_reads, gene_body_reads, pausing_index))
    return(
        list(refseq_tss_coord, tss_reads, gene_body_reads, pausing_index))
}

calcAntiSensePauseRatio2 <- compiler::cmpfun(calcAntiSensePauseRatio)

## lineprof(
calcAntiSensePauseRatio2(data = neg_dat,
                         refseq_tss_coord = 39549838,
                         tss_range = 1000,
                         tss_win_size = 50,
                tss_calc_step = 5,
                body_threshold = 5,
                body_cutoff = 2500)
## )

## ## TODO Save File, Visualize...
## print("[POS STRAND] Executing in parallel")
## pausing_indices <- mclapply(tss_test$start, function(x)
##     calcSensePauseRatio2(data = pos_dat,
##                          refseq_tss_coord = x,
##                     tss_range = 1000,
##                     tss_win_size = 50,
##                     tss_calc_step = 5,
##                     body_threshold = 5,
##                     body_cutoff = 2500),
##     mc.cores = no_cores
##     )

## print("[NEG STRAND] Executing in parallel")
## pausing_indices <- mclapply(tss_test$start, function(x)
##     calcAntiSensePauseRatio2(data = neg_dat,
##                              refseq_tss_coord = x,
##                              tss_range = 1000,
##                              tss_win_size = 50,
##                              tss_calc_step = 5,
##                              body_threshold = 5,
##                              body_cutoff = 2500),
##     mc.cores = no_cores
##     )

message("Executing linearly")
pausing_indices <- lapply(tss_test$start, function(x)
    calcAntiSensePauseRatio2(data = neg_dat,
                             refseq_tss_coord = x,
                             tss_range = 1000,
                         tss_win_size = 50,
                         tss_calc_step = 5,
                    body_threshold = 5,
                    body_cutoff = 2500)
    )

message("Finished parsing pausing indices!")

save(pausing_indices, file="/scratch/Users/zama8258/PAUSING_PRELIM.Rda")

## Generate Figures
## f <- read_tsv("/scratch/Users/zama8258/PAUSING_PRELIM.tsv", col_names = FALSE)
f <- load("/scratch/Users/zama8258/PAUSING_PRELIM.Rda")

dd <- as.tibble(do.call(rbind, pausing_indices))
dd$V1 <- unlist(dd$V1)
dd$V2 <- unlist(dd$V2)
dd$V3 <- unlist(dd$V3)
dd$V4 <- unlist(dd$V4)

## ddf <- filter(dd, !is.na(dd$V4))

ddf <- filter(dd, V4 > 10)

library(ggthemes)
library(ggsci)
library(wesanderson)

ggplot(data = ddf, mapping=aes(x=V4), alpha=1/10) +
    geom_histogram(binwidth=5) + theme_tufte() +
    labs(x = "Pausing Index",
         y= "Count",
         title= "Distribution of Pausing Indices") 
ggsave("pause_plot_neg_prelim.png", plot = last_plot(), device = "png")

## Clean up the cluster once we're done.
## stopCluster(cl)
