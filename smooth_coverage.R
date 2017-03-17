suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
library(mmand)
library(assertthat)

# Given a file of read coverage, smooth it before
# calling the no tetrad info HMM.
# Must specify a two args (the coverage file and out file)
if (!exists("argv")) argv <- commandArgs(trailingOnly=TRUE)
assert_that(length(argv) == 2)
filename <- argv[1]
outfile <- argv[2]

cols <- c('chrom', 'start', 'end', 'bin_name', 'y')
dat <- fread(filename, col.names=cols)
sizes <- dat$end - dat$start
assert_that(all(sizes > 0))
binsize <- median(sizes)

# median filter
n <- 11
k <- shapeKernel(n)
dat.filtered <- group_by(dat, chrom) %>%
    mutate(medfilt=round(medianFilter(y, k)))
write.table(select(dat.filtered, -y), file=outfile,
    quote=F, row.names=F, col.names=F, sep="\t")
