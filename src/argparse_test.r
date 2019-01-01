suppressMessages(library("tidyverse"))
## Parse arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if (length(args) < 1) {
    args <- c("--help")
}

## Help section
if ("--help" %in% args) {
    cat("
      R Argument Parsing
 
      Arguments:
      --srr=SRR000000    - SRR Number
      --help              - print this text
 
      Example:
      ./script.r --srr=SRR123456 \n\n")
    
    q(save = "no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))

if (is.null(argsL)) {
    print("You must provide an SRR")
    q()
}

## Make sure we parse the correct srr
message(str_c("SRR is: ", argsL))

