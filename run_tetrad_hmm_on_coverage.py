#!/usr/bin/env python

"""Take BED files of coverage across four spores.
Read into a list of numpy matrices, one numpy matrix per chromosome.

Run a ContinuousMultinomialHMM hmmlearn model on this data.

This version estimates emission probabilities rather than
getting them beforehand via EM estimation.

Pass around random seed correctly as explained by Robert Kern in
https://stackoverflow.com/questions/5836335/consistenly-create-same-random-numpy-array
"""
import sys, argparse, os, re, copy, itertools, pickle
import numpy as np
from numpy.random import RandomState
import pandas as pd
from hmmlearn import hmm
from continuous_multinomial_hmm import ContinuousMultinomialHMM
from scipy.stats import describe

def main():
    args = getArgs()
    prng = RandomState(args.seed)

    # First, read in BED files of coverage in all four spores.
    # Next find the bin size. Should be the same for all bins and all BED files
    cov = readFourCoverageBed(args.bedfiles)
    # cov is a pandas groupby object. Can say, e.g. cov.get_group('chrI').head()
    index1 = cov.groups.values()[0][0]
    windowsize = index1[2] - index1[1]
    for chrom, df in cov:
        index = df.index
        sizes = index.get_level_values('end') - index.get_level_values('start')
        assert np.all(sizes == windowsize)

    # calculate recombination fraction between windows
    # c == recomb rate per base (.01/2000) * nbases
    if args.verbose:
        print >> sys.stderr, "calculated c using window size %dbp" % windowsize
    c = 5e-06*windowsize

    X, indices = cov2numpy(cov, args)
    # indices is a list of len == n chroms. Each chrom item is a list of
    # tuples that are (chrom, binStart, binEnd)
    lengths = [xx.shape[0] for xx in X]
    Xcat = np.concatenate(X)

    # Initialize the HMM, then run it on all data.
    # Run on separate, randomly initialized, emission prob matrices
    model = initialize_hmm(c, prng=prng, LOH=args.include_LOH_state)
    model_runs = run_hmm(model, Xcat, lengths, args=args, prng=prng)
    if args.verbose:
        print >> sys.stderr, "Using %d good HMM runs" % len(model_runs)
    if args.save_emission_matrices_file is not None:
        matrices = dict(('matrix%d' % i, run.emissionprob_) \
            for i, run in enumerate(model_runs))
        np.savez(args.save_emission_matrices_file, **matrices)
    best_run = checkRunResults(model_runs, Xcat, lengths)

    # use predict_proba() to get most likely ancestry for each window,
    #     on each chrom, for this tetrad.
    if args.verbose:
        print >> sys.stderr, "Calculating posterior probs for each state"
    probs = best_run.predict_proba(Xcat, lengths=lengths)

    # get state order for emission probability matrix.
    # usually it is 1100, 0101, 0110, 1001, 1010, 0011
    # because I start it near there when initializing.
    # But it doesn't have to be in that order.
    roundmat = np.empty(shape=best_run.emissionprob_.shape, dtype=np.int)
    for i, row in enumerate(best_run.emissionprob_):
        ranks = row.argsort().argsort()
        roundmat[i, :] = ranks
    roundmat[roundmat < 2] = 0
    roundmat[roundmat > 1] = 1
    # should be 2 across all rows:
    assert np.allclose(roundmat.sum(axis=1), \
        np.full((roundmat.shape[0], ), 2, dtype=np.int))
    s2r = dict()    # state to row dict,
                    # key=state eg.1001; val=row in emission probability matrix
    for i in range(6):  # only look at the estimated components even if LOH state included
        state = ''.join(["%d" % j for j in roundmat[i, :]])
        assert state not in s2r
        s2r[state] = i

    # row[i] is posterior probabilities for each state of bin $i$
    # (L->R in row is L->R/top->bottom on trprobs matrix)
    #
    # add together the probabilities of the states representing coverage
    #     of that strain (3/6 states) to get a prob of strain ancestry
    #     for each window. Should sum to 2 for any specific window, across
    #     all four spores.
    # normal state order is:
    # state 1100, 0101, 0110, 1001, 1010, 0011, (+LOH)
    # which would mean:
    # spore 1, ancestry from this strain, is sum of probs[, [0, 3, 4]]
    # spore 2, ancestry from this strain, is sum of probs[, [0, 1, 2]]
    # spore 3, ancestry from this strain, is sum of probs[, [2, 4, 5]]
    # spore 4, ancestry from this strain, is sum of probs[, [1, 3, 5]]
    onerows = [s2r['1100'], s2r['1001'], s2r['1010']]
    tworows = [s2r['1100'], s2r['0101'], s2r['0110']]
    threerows = [s2r['0110'], s2r['1010'], s2r['0011']]
    fourrows = [s2r['0101'], s2r['1001'], s2r['0011']]
    ancestry = pd.DataFrame({'one':probs[:, onerows].sum(axis=1),
        'two'   : probs[:, tworows].sum(axis=1),
        'three' : probs[:, threerows].sum(axis=1),
        'four'  : probs[:, fourrows].sum(axis=1)},
        columns=['one', 'two', 'three', 'four'],
        index=pd.MultiIndex.from_tuples(np.concatenate(indices)))
    if not args.include_LOH_state:
        # should be 2 across all rows if no LOH state:
        assert np.allclose(ancestry.sum(axis=1), \
            np.full((ancestry.shape[0],), 2, dtype=float))
    else:
        # LOH means ancestry from the non-ref strain across all spores
        lohprob = probs[:, 6]
        for key in ancestry:
            ancestry[key] += lohprob
        loh_highprob = lohprob > 0.5
        assert np.allclose(ancestry.ix[~loh_highprob].sum(axis=1), \
            np.full(((~loh_highprob).sum(),), 2, dtype=float))
        assert np.allclose(ancestry.ix[loh_highprob].sum(axis=1), \
            np.full((loh_highprob.sum(), ), 4, dtype=float))

    # Print a little table that tells us, for each state,
    # how many bins/windows are in that state with a posterior probability
    # >0.9.
    probs = pd.DataFrame(probs)
    maxprobs = probs.max(axis=1) > 0.9
    idx = probs.ix[maxprobs, :].idxmax(axis=1)   # index of prob > 0.9
    if args.verbose:
        print >> sys.stderr, "Number of bins of high prob for each state:"
        print >> sys.stderr, idx.value_counts().sort_index()

    np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
    if not args.hide_emissionmat:
        print "# emission probability matrix:"
        for row in model.emissionprob_:
            print "# %s" % row
    if not args.hide_transitionmat:
        print "# transition probability matrix:"
        for row in model.transmat_:
            print "# %s" % row
    if args.include_LOH_state:
        print "# LOH probability is included in each of the four spores %s" \
            % ("(assumes LOH for non-reference sequence)")

    print ancestry.to_csv(None, sep="\t", header=False, columns=None, \
        index=True, float_format="%.4g"),

def anomalous_runlen_filtering(bools, n):
    """Given a boolean pd.Series marking loci with 'anomalous' data
    (i.e. doesn't fit any pattern we expect, maybe due to repetitive
    sequence, etc), convert it into a new series with True only
    at the indices with <= n consecutive Trues.
    """
    bb = bools.copy()
    newbool = np.zeros(bb.shape, dtype=np.bool)
    tf_counts = [(b, len(list(g))) for b, g in itertools.groupby(bb)]
    lens = np.array([tup[1] for tup in tf_counts])
    csum = np.cumsum(lens)
    starts = np.concatenate([np.array([0]), csum[:-1]])
    ends = csum
    t_indices = np.array([tup[0] for tup in tf_counts])
    t_runs_short = (lens <= n) & t_indices
    for tup in zip(starts[t_runs_short], ends[t_runs_short]):
        newbool[tup[0]:tup[1]] = True
    return newbool

def checkRunResults(model_runs, Xcat, lengths, atol=0.005, rtol=0.03):
    """This function checks HMM run results.
    I have a list of HMM runs where all runs had fairly good loglik.
    Now I need to check that the inferred transition and emission probability
    matrices are reasonably close for all runs.

    Return best run (as measured by loglik). If matrices are not reasonably
    close, raise an error.
    """
    logliks = [m.score(Xcat, lengths=lengths) for m in model_runs]
    best_run = np.argmax(logliks)
    best_trprobs = model_runs[best_run].transmat_
    best_eprobs = model_runs[best_run].emissionprob_

    for mod in model_runs:
        assert np.allclose(best_trprobs, mod.transmat_, atol=atol, rtol=rtol)
        assert np.allclose(best_eprobs, mod.emissionprob_, atol=atol, rtol=rtol)
    return model_runs[best_run]

def cov2numpy(cov, args, verbose=True):
    """Convert the pandas datasets to numpy.
    Ignore sites subject to args.anomalous_threshold"""
    if args.anomalous_threshold < 1 and args.anomalous_threshold != 0:
        is_anomalous = three_spores_fraction_above
    else:
        is_anomalous = three_spores_readcount_above

    X, Xnames, indices = list(), list(), list()
    anomalous_sites, all_sites = 0, 0
    for chrom, dataframe in cov:
        Xnames.append(chrom)
        dat = dataframe.copy()
        indices.append(dat.index)

        anomalous_sites_bool = dat.apply(is_anomalous, axis=1, \
            n=args.anomalous_threshold)
        anomalous_sites_bool_runlen = anomalous_runlen_filtering( \
            anomalous_sites_bool, args.anomalous_runlen)
        anomalous_sites += anomalous_sites_bool_runlen.sum()
        all_sites += len(anomalous_sites_bool_runlen)
        dat.loc[anomalous_sites_bool_runlen, :] = 0

        X.append(np.matrix(dat))
    if verbose:
        print >> sys.stderr, "Filtered out %d anomalous sites" % anomalous_sites
        print >> sys.stderr, "Out of %d total sites" % all_sites
    return X, indices

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("bedfiles", nargs=4, metavar="file", \
        help="List of BED files of coverage, one for each of " + \
        "the four spores. Four columns, tab-delimited: " + \
        "chrom, start, end, coverage")
    parser.add_argument("--anomalous_runlen", default=2, type=int, \
        help="Do not filter out anomalously thresholded sites " + \
        "if there are more than this number consecutively occurring." + \
        " (useful to set high if there is LOH data)", metavar='N')
    parser.add_argument("--anomalous_threshold", default=0.15, type=float, \
        help="Number of reads or min fraction of reads present in >=3" + \
        " spores above which the window should be filtered out", metavar="N")
    parser.add_argument("--hide_emissionmat", action="store_true", \
        help="Do not print emission matrix at top of output file")
    parser.add_argument("--hide_transitionmat", action="store_true", \
        help="Do not print transition matrix at top of output file")
    parser.add_argument("--include_LOH_state", action="store_true", \
        default=False, help="include state with expr from only one ancestor")
    parser.add_argument("--min_good_reps", type=int, default=None, metavar="N",\
        help="number of reps with reasonably good loglik required")
    parser.add_argument("--nreps", type=int, default=12, \
        help="number of times the HMM is fit with different emission" + \
        " probability matrix random initializations", metavar="N")
    parser.add_argument("--save_emission_matrices_file", default=None, \
        metavar="emat.npz", \
        help="Name of file to save emission matrices in np.savez format")
    parser.add_argument("--seed", help="seed for HMM", type=int, \
        default=1234, metavar="N")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    if args.min_good_reps is None:
        args.min_good_reps = args.nreps - 4
    assert args.min_good_reps >= 3
    assert args.nreps >= args.min_good_reps
    if args.save_emission_matrices_file is not None:
        assert args.save_emission_matrices_file.endswith('.npz')
    return args

def getBedLine(filehandle):
    line = filehandle.readline().rstrip('\n')
    if line == '': return ''
    words = line.split('\t')
    chrom, start, end = words[0], int(words[1]), int(words[2])
    cov = float(words[-1])
    return (chrom, start, end, cov)

def pandasReadBedCoverage(filename, colname='coverage'):
    dat1 = pd.read_table(filename, comment='#', header=None,
        names=['chrom', 'start', 'end', colname])
    # make a MultiIndex for chrom:start-end.
    # False == do not return index, None == do not name tuples
    index = [t for t in dat1[['chrom', 'start', 'end']].itertuples(False, None)]
    dat1.drop(['chrom', 'start', 'end'], axis=1, inplace=True)
    dat1.index = pd.MultiIndex.from_tuples(index,
        names=['chrom', 'start', 'end'])
    return dat1

def initialize_hmm(c, prng, LOH, n_iter=1000, tol=0.01):
    """Initialize the HMM.
    c = recombination rate per window size. E.g. 0.01 per 2kb (1cM/2kb),
        or 5e-04 per 100bp.
    components = hidden states
    features = symbols emitted by each state

    Based in part on complex_tetrad_example.py ...
    """
    # transition matrix. shape (n, n)
    nc = 7 if LOH else 6
    s = 1 - c*4     # "s"elf == no recombination
                    # four possible ways to switch states by recombination
                    # one "inaccessible" state would require two recomb. events
                    # 12/1/16 add an LOH state (bottom row & right row)
                    # only enter LOH at beginning of chrom, then stuck there
    if c >= 0.2:
        # the above scheme fails otherwise since it leads to some negative probs
        # if c is high we want to the prob of transitioning to all other
        # states to be equivalent. We probably won't much trust results with
        # bins this large, since they are likely to contain multiple recombination
        # breakpoints.
        s, c = 0.2, 0.2
                      # to state (read bottom to top)
                      #  0  1  0  1  0  1  1
                      #  0  0  1  0  1  1  1
                      #  1  1  1  0  0  0  1
                      #  1  0  0  1  1  0  1
    trprobs = np.array([[s, c, c, c, c, 0, 0],  # from state 1100
                        [c, s, c, c, 0, c, 0],  # from state 0101
                        [c, c, s, 0, c, c, 0],  # from state 0110
                        [c, c, 0, s, c, c, 0],  # from state 1001
                        [c, 0, c, c, s, c, 0],  # from state 1010
                        [0, c, c, c, c, s, 0],  # from state 0011
                        [0, 0, 0, 0, 0, 0, 1]]) # from state 1111
    if not LOH:
        inds = np.arange(6)
        trprobs = trprobs[inds[:, None], inds]
    start = np.repeat(np.array([1./nc]), nc)			# shape (nc,)

    # params and init_params can each be 'ste' or any subset.
    # init_params determines which parameters should be initialized
    # automatically, see http://hmmlearn.readthedocs.io/en/latest/api.html
    # and https://github.com/hmmlearn/hmmlearn/issues/62
    # params is which parameters should be estimated when fit() is called
    # i.e. updated during the training process
    # 's'=startprob, 't'=transition probs, 'e'=emission probs
    model = ContinuousMultinomialHMM(n_components=nc, random_state=prng,
        n_iter=n_iter, tol=tol, verbose=True, init_params='', params='te',
        est_transition_indices=np.arange(6))
    model.transmat_ = trprobs
    model.startprob_ = start
    return model

def initialize_random_emission_probs(prng, LOH, error_rate=0.002):
    """Initialize an emission probability matrix with random probabilities.
    Keep it subject to the structure noted below, where the two probabilities
    that are highest in each row are in fixed locations in the matrix.
    Here is an example deterministic matrix initialization:
    error_rate = 0.002
    e = error_rate/2
    t = (1 - error_rate)/2  # "t"rue data means no mismapping/seq error
                            # emission probs shape (n_components, n_features)
                            # emit state (read bottom to top)
    if LOH, there are two even states.
    even = [0.25, 0.25, 0.25, 0.25]
    emissionp_mat = np.array([[t, t, e, e],     # in state 1100
                              [e, t, e, t],     # in state 0101
                              [e, t, t, e],     # in state 0110
                              [t, e, e, t],     # in state 1001
                              [t, e, t, e],     # in state 1010
                              [e, e, t, t],     # in state 0011
                              even])            # in state 1111
    """
    # record which columns are the two higher probs and which are the two lower
    nrows = 7 if LOH else 6
    mat = np.empty((nrows, 4))
    # highprobs = [[0, 1], [1, 3], [1, 2], [0, 3], [0, 2], [2, 3]]
    # lowprobs = [[2, 3], [0, 2], [0, 3], [1, 2], [1, 3], [0, 1]]
    # for i in range(6):
    #     fournums = prng.random_sample(4)
    #     fournums = sorted(fournums/sum(fournums), reverse=True) # biggest first
    #     mat[i, highprobs[i]] = prng.choice(fournums[:2], size=2, replace=False)
    #     mat[i, lowprobs[i]] = prng.choice(fournums[2:], size=2, replace=False)
    mat[0, ] = np.random.dirichlet(alpha=[10, 10, 1, 1])
    mat[1, ] = np.random.dirichlet(alpha=[1, 10, 1, 10])
    mat[2, ] = np.random.dirichlet(alpha=[1, 10, 10, 1])
    mat[3, ] = np.random.dirichlet(alpha=[10, 1, 1, 10])
    mat[4, ] = np.random.dirichlet(alpha=[10, 1, 10, 1])
    mat[5, ] = np.random.dirichlet(alpha=[1, 1, 10, 10])
    if LOH:
        even = [0.25, 0.25, 0.25, 0.25]
        mat[6, ] = even
    assert np.allclose(mat.sum(axis=1), 1)
    return mat

def readFourCoverageBed(filenames):
    "Return in a bunch of coverage data in BED format, 1 for each of 4 spores"
    assert len(filenames) == 4
    datlist = [pandasReadBedCoverage(filenames[i], i+1) for i in range(4)]
    # align all datasets using indexes.
    alldat = pd.concat(datlist, axis=1)
    # assert no missing data
    assert not np.any(pd.isnull(alldat))
    # return data split by chrom
    # alldat.xs(chrom, level=0) would give data for one chrom
    # alldat.xs(start, level=1) would give data for all bins starting@start
    # alldat.xs((chrom, end), level=(0, 2)) works for finer interrogation
    # can also use alldat.loc[] but that may require lexsorting:
    # alldat.sort_index(inplace=True)
    # alldat.loc[(slice(None), 100, slice(None)), :]
    return alldat.groupby(level='chrom', sort=False)

def readFourLines(handles):
    assert len(handles) == 4
    lines = list()
    for i in range(len(handles)):
        lines.append(handles[i].readline().rstrip('\n'))
    if lines[0] == '':
        assert lines[1] == lines[2] == lines[3] == ''
        return None
    return lines

def run_hmm(model, Xcat, lengths, args, prng, verbose=True):
    """Initialize the HMM, then run it on all data.
    Run on separate, randomly initialized, emission prob matrices

    1) run HMM args.nreps times
    2) obtain at least arg.min_good_reps that have reasonably good
       likelihood. Our goal is to compare the results of these good
       runs and make sure their transition and emission probability
       matrices are reasonably similar ... i.e. the results are stable.
       How do we define good runs? The goal is to avoid including really
       bad runs, but if very good runs are really rare, we might never
       get enough really good runs. So we find the difference between the
       best and third-best runs (logdiff), and we take all runs within
       (logdiff + eps) of the best run.

    Return a list of the fitted models.
    """
    orig_trprobs = model.transmat_.copy()
    model_list, logliks = list(), list()
    n = 1
    for i in range(args.nreps):
        if verbose: print >> sys.stderr, "HMM run number %d" % n
        n += 1
        emissionp_mat = initialize_random_emission_probs(prng, LOH=args.include_LOH_state)
        model.emissionprob_ = emissionp_mat
        model.transmat_ = orig_trprobs
        model = model.fit(Xcat, lengths=lengths)
        assert model.monitor_.converged
        model_list.append(model)
        logliks.append(model.score(Xcat, lengths=lengths))

    best_run = np.argmax(logliks)
    loglik_best = logliks[best_run]
    # set eps to be within 5% of logdiff
    eps = np.log(1.05)  # exp(loglik_best)/exp(this_loglik) < eps means this_loglik is almost as good, or better
    logdiff = loglik_best - sorted(logliks, reverse=True)[2]
    models, model_logliks = list(), list()
    for i in range(args.nreps):
        if loglik_best - logliks[i] < logdiff + eps:
            models.append(model_list[i])
            model_logliks.append(logliks[i])

    while len(models) < args.min_good_reps:
        if verbose: print >> sys.stderr, "HMM run number %d" % n
        n += 1
        if n > 500:
            print >> sys.stderr, "ERROR: ran HMM 500 times without achieving satisfactory convergence. What is wrong?"
            print >> sys.stderr, "logliks for first %d reps were %s" % (args.nreps, ', '.join(["%.3g" % l for l in sorted(model_logliks)]))
            print >> sys.stderr, "loglik for last model run was %.3g" % this_loglik
            sys.exit(1)
        emissionp_mat = initialize_random_emission_probs(prng, LOH=args.include_LOH_state)
        model.emissionprob_ = emissionp_mat
        model.transmat_ = orig_trprobs
        model = model.fit(Xcat, lengths=lengths)
        this_loglik = model.score(Xcat, lengths=lengths)
        # keep only if reasonably close in loglik:
        if abs(this_loglik - loglik_best) < logdiff + eps:
            models.append(model)
    return models

def three_spores_fraction_above(row, n, min_reads=2):
    "return True if 3 or more spores has at least fraction n of the total reads"
    assert 0 < n < 1
    if all(row <= min_reads): return False
    return sum(row/sum(row) > n) > 2

def three_spores_readcount_above(row, n):
    "return True if 3 or more spores have read count above n"
    return sum(row > n) > 2

if __name__ == "__main__":
    sys.exit(main())
