#!/usr/bin/env Rscript --vanilla 

suppressPackageStartupMessages(require(optparse)) 

option_list = list(
  make_option(c("-p", "--pval"), 
              type = "character", 
              default=NULL,
              help = "Parent ppi file"),
  make_option(c("-o", "--out_pdf"), 
              type = "character", 
              default=NULL,
              help = "Name of pdf output")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

renv::load(here::here())

pvals <- read.table(here::here(opt$pval), header=T)

observed <- sort(pvals$GenePvalue)
lobs <- -(log10(observed))

expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))


pdf(here::here(opt$out_pdf), width=6, height=6)
plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()


