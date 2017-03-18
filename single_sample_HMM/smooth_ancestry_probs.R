suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
library(mmand)
library(assertthat)

# Given a file of ancestry probabilities (i.e. output of single-sample HMM),
# smooth the probabilities to 0/1 and filter out small ancestry bins
# Must specify a two args (the coverage file and out file)
if (!exists("argv")) argv <- commandArgs(trailingOnly=TRUE)
assert_that(length(argv) == 2)
filename <- argv[1]
outfile <- argv[2]

cols <- c('chrom', 'start', 'end', 'ref', 'nonref')
dat <- fread(filename, col.names=cols) %>%
  mutate(round_ref=round(ref), round_nonref=round(nonref),
    sum_rprobs=round_ref + round_nonref)
nonref <- select(dat, -ref, -nonref, -round_ref) %>%
  rename(nonref=round_nonref)

namebins <- function(hardcutoff) {
  # function from no_tetrad_info/call_breakpoints.R
  r <- rle(hardcutoff)
  rep(1:length(r$lengths), r$lengths)
}
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

min_bin_size <- 20001
a <- group_by(nonref, chrom) %>% mutate(bin_name=namebins(nonref)) %>%
  group_by(chrom, bin_name) %>% mutate(bin_size=max(end) - min(start)) %>%
  mutate(too_small=bin_size < min_bin_size | sum_rprobs != 1, newcall=nonref) %>%
  arrange(chrom, start, end)
a$newcall[a$too_small] <- -17
a2 <- group_by(a, chrom) %>% mutate(final_call=fix_missing(newcall)) %>%
  ungroup() %>% arrange(chrom, start, end) %>%
  select(chrom, start, end, final_call)

write.table(a2, file=outfile, quote=F, row.names=F, col.names=F, sep="\t")
