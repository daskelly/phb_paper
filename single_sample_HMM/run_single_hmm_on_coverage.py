#!/usr/bin/env python

"""Take a BED file of coverage across one strain,
and BED files that give read mapping probabilities
(obtained from simulation) for each bin.

Read into a list of numpy vectors, one numpy vector per chromosome.

Run a ContinuousBinomialHMM hmmlearn model on this data.
Use the predict_proba() method to use the forward-backward
algorithm to obtain probabilities of ancestry for each bin.

Pass around random seed correctly as explained by Robert Kern in
https://stackoverflow.com/questions/5836335/consistenly-create-same-random-numpy-array
"""
import sys, argparse, os, re, copy, itertools, pickle
import numpy as np
from numpy.random import RandomState
import pandas as pd
from hmmlearn import hmm
from continuous_binomial_hmm import ContinuousBinomialHMM

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bedfile", metavar="coverage.bed", \
        help="BED file of coverage. Five columns, tab-delimited: " + \
        "chrom, start, end, name, coverage")
    parser.add_argument("sim_p_ancestry", metavar="p_ancestry.bed", \
        help="BED file: chrom, start, end, name, cov, frac")
    parser.add_argument("sim_p_no_ancestry", metavar="p_no_ancestry.bed", \
        help="same format as sim_p_ancestry")
    parser.add_argument("--seed", help="seed for HMM", type=int, \
        default=1234, metavar="N")
    parser.add_argument("--summary_state_minprob", default=0.9, type=float, \
        metavar="f")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    prng = RandomState(args.seed)

    # First, read in BED file of coverage.
    coverage = pandasReadBedCoverage(args.bedfile, extracols=['name', 'coverage'])
    coverage.drop('name', inplace=True, axis=1)
    sizes = coverage.index.get_level_values('end') - coverage.index.get_level_values('start')
    windowsize = int(np.round(np.median(sizes)))
    if args.verbose:
        print >> sys.stderr, "Using a window size of {:d}bp (median bin size)".format(windowsize)

    # Next read in bin-specific probabilities
    # these are obtained via simulation
    sim_p_ancestry = pandasReadBedCoverage(args.sim_p_ancestry, \
        extracols=['name', 'cov', 'p_ancestry'])
    sim_p_ancestry.drop(['name', 'cov'], inplace=True, axis=1)
    sim_p_no_ancestry = pandasReadBedCoverage(args.sim_p_no_ancestry, \
        extracols=['name', 'cov', 'p_no_ancestry'])
    sim_p_no_ancestry.drop(['name', 'cov'], inplace=True, axis=1)
    # merge dataframes to be sure bins are aligned properly...
    alldat = pd.concat([coverage, sim_p_ancestry, sim_p_no_ancestry], axis=1)
    alldat.sort_index(level=['chrom', 'start', 'end'], inplace=True)

    # calculate recombination fraction between windows
    # c == recomb rate per base (.01/2000) * nbases
    if args.verbose:
        print >> sys.stderr, "calculated c using window size {:d}bp".format(windowsize)
    c = 5e-06*windowsize

    cov = alldat.groupby(level='chrom', sort=False)
    X, indices = cov2numpy(cov)
    # indices is a list of len == n chroms. Each chrom item is a list of
    # tuples that are (chrom, binStart, binEnd)
    lengths = [xx.shape[0] for xx in X]
    Xcat = np.concatenate(X)
    N = np.repeat(alldat['coverage'].sum(), alldat.shape[0])
    P = np.array((alldat['p_no_ancestry'], alldat['p_ancestry']))

    # Initialize the HMM, then compute posterior probabilities
    model = initialize_hmm(nbins=alldat.shape[0], c=c, prng=prng)
    Xcat = np.vstack((Xcat, N, P)).T
    probs = model.predict_proba(Xcat, lengths=lengths)

    # Print a little table that tells us, for each state,
    # how many bins/windows are in that state with a posterior probability
    # >0.9.
    probs = pd.DataFrame(probs)
    probs.columns = ['ref', 'non-ref']
    probs.index = alldat.index
    maxprobs = probs.max(axis=1) > args.summary_state_minprob
    idx = probs.ix[maxprobs, :].idxmax(axis=1)   # index of prob > 0.9
    if args.verbose:
        print >> sys.stderr, "Looking at {:d} total bins".format(alldat.shape[0])
        print >> sys.stderr, "Number of bins of high prob for each state:"
        print >> sys.stderr, idx.value_counts().sort_index()

    print "# ",
    print probs.to_csv(None, sep="\t", header=True, columns=None, \
        index=True, float_format="%.4g"),

def cov2numpy(cov, verbose=True):
    """Convert the pandas datasets to numpy."""
    X, indices = list(), list()
    for chrom, dataframe in cov:
        X.append(np.array(dataframe['coverage']))
        indices.append(dataframe.index)
    return X, indices

def pandasReadBedCoverage(filename, extracols=['coverage']):
    dat1 = pd.read_table(filename, comment='#', header=None,
        names=['chrom', 'start', 'end'] + extracols)
    # make a MultiIndex for chrom:start-end.
    # False == do not return index, None == do not name tuples
    index = [t for t in dat1[['chrom', 'start', 'end']].itertuples(False, None)]
    dat1.drop(['chrom', 'start', 'end'], axis=1, inplace=True)
    dat1.index = pd.MultiIndex.from_tuples(index,
        names=['chrom', 'start', 'end'])
    return dat1

def initialize_hmm(nbins, c, prng, n_iter=1000, tol=0.01):
    """Initialize the HMM.
    c = recombination rate per window size. E.g. 0.01 per 2kb (1cM/2kb),
        or 5e-04 per 100bp.
    components = hidden states
    features = symbols emitted by each state (1 for binomial data since 1-dimensional)
    """
    s = 1 - c     # "s"elf == no recombination
    # states (L-R or top-bottom) are [ref, strain]
    trprobs = np.array([[s, c],
                        [c, s]])
    nc = 2  # number of components == hidden underlying states == [ref, strain]
    start = np.repeat(np.array([1./nc]), nc)			# shape (nc,)

    # emission probs matrix is not used. This is because we
    # pass emission probs along with the data X to predict_proba().
    # The reason is because we have emission probs that vary by bin
    # and need to be matched up with data points in X properly.
    # See _compute_log_likelihood() in continuous_binomial_hmm.py
    emissionp_mat = np.zeros((nbins, nc))

    # params and init_params can each be 'ste' or any subset.
    # init_params determines which parameters should be initialized
    # automatically, see http://hmmlearn.readthedocs.io/en/latest/api.html
    # params is which parameters should be estimated when fit() is called
    # 's'=startprob, 't'=transition probs, 'e'=emission probs
    model = ContinuousBinomialHMM(n_components=nc, random_state=prng,
        n_iter=n_iter, tol=tol, verbose=True, init_params='', params='')
    model.transmat_ = trprobs
    model.startprob_ = start
    model.emissionprob_ = emissionp_mat
    return model

if __name__ == "__main__":
    sys.exit(main())
