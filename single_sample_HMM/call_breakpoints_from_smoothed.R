suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
library(mmand)
library(assertthat)
library(tidyr)
suppressPackageStartupMessages(library(tibble))

# Must specify two args. First arg is the output of
# the "no tetrad info" version of run_hmm_on_coverage.py BUT
# this file has been run through smooth_ancestry_probs.R.
# The file has format:
# chrom start end nonref_prob
# Second arg is the name of the output file

if (!exists("argv")) argv <- commandArgs(trailingOnly=TRUE)
assert_that(length(argv) == 2)
outfile <- argv[2]
filename <- argv[1]
assert_that(endsWith(filename, ".bed"))
bfile <- basename(filename)
prefix <- substring(bfile, 1, nchar(bfile)-4)
cols <- c('chrom', 'start', 'end', 'nonref_prob')
dat <- as_tibble(fread(filename, col.names=cols))

fix_missing <- function(vals) {
  # We have a vector of 0/1/-17.
  # Any -17 between 0 and 0 are set to 0
  # Any -17 between 1 and 1 are set to 1
  # Any -17 between 0 and 1, we set a gradual range in between
  # Any -17 at start or end of chrom, we set to the first non -17 value.
  newvals <- vals
  r <- rle(vals)
  n <- length(r$lengths)
  idx_ends <- cumsum(r$lengths)
  idx_starts <- c(1, idx_ends[-length(idx_ends)]+1)
  if (r$values[1] == -17) {
    newvals[idx_starts[1]:idx_ends[1]] <- r$values[2]
    r$values[1] <- r$values[2]
  }
  if (r$values[n] == -17) {
    newvals[idx_starts[n]:idx_ends[n]] <- r$values[n-1]
    r$values[n] <- r$values[n-1]
  }
  for (i in 1:n) {
    if (r$values[i] == -17) {
      val_before <- r$values[i-1]
      val_after <- r$values[i+1]
      if (val_before == val_after) {
        newvals[idx_starts[i]:idx_ends[i]] <- val_before
      } else {
        nn <- length(idx_starts[i]:idx_ends[i])
        ss <- seq(val_before, val_after, length.out=nn+2)[2:(nn+1)]
        newvals[idx_starts[i]:idx_ends[i]] <- ss
      }
    }
  }
  assert_that(all(newvals >= 0 & newvals <= 1))
  newvals
}

find.transitions <- function(df) {
  x <- df$nonref_prob
  x[x != 0 & x != 1] <- -17
  r <- rle(x)
  n <- length(r$lengths)
  if (n == 1) return(NULL)
  idx_ends <- cumsum(r$lengths)
  idx_starts <- c(1, idx_ends[-length(idx_ends)]+1)

  assert_that(r$values[1] != -17)
  assert_that(r$values[n] != -17)
  # need to do the following
  # for 0/1 or 1/0 transitions, just take the end of the 'before' bin and
  # start of the 'after' bin and average
  # for 0/-17/1 or 1/-17/0, take the end of the first bin and start of
  # the last bin and average.
  n_transitions <- sum(r$values != -17) - 1
  starts <- numeric(n_transitions)
  ends <- numeric(n_transitions)
  j <- 1
  for (i in 1:(n-1)) {
    if (r$values[i] == -17) next
    starts[j] <- df$end[idx_ends[i]]
    if (j > 1) ends[j-1] <- df$start[idx_starts[i]]
    j <- j + 1
  }
  ends[j-1] <- df$start[idx_starts[n]]
  newdat <- tibble(chrom=df$chrom[1], pos=rowMeans(cbind(starts, ends)))
  return(newdat)
}

# nonref_prob is already smoothed, so I just need to find 0/1 and 1/0
# transitions, and put the transition in between.
d2 <- group_by(dat, chrom) %>% dplyr::do(as.data.frame(find.transitions(.))) %>%
  mutate(end=round(pos), start=end - 1) %>% select(-pos) %>%
  select(chrom, start, end)

fwrite(d2, sep="\t", col.names=FALSE, file=outfile)
