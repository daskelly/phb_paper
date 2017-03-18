"""Extend hmmlearn with an HMM that emits binomial distributed
data that is conditional on total coverage across the genome
and has bin-specific probabilities for ancestry/no ancestry."""
import numpy as np
import sys
from hmmlearn import base
from scipy.stats import binom

class ContinuousBinomialHMM(base._BaseHMM):
    """Extends base._BaseHMM in the hmmlearn package for use with data
    where the emissions are binomial data conditional on a given
    sample size at each site (sequence coverage)."""
    def __init__(self, n_components=2, algorithm="viterbi", random_state=None,
        n_iter=1000, tol=1e-2, verbose=False, params="ste", init_params="ste"):
        base._BaseHMM.__init__(self, n_components, algorithm=algorithm,
        random_state=random_state, n_iter=n_iter, tol=tol, verbose=verbose,
        params=params, init_params=init_params)

    def _compute_log_likelihood(self, X):
        """Return the log of the binomial probability density.
        Needs to return a matrix that is of shape
        (n_obs in sequence, n_components).
        In order to accommodate having variable probabilities in each bin
        along the genome (and having the (X, lengths) model of hmmlearn),
        I need to synthetically combine X and p so that they are divided
        up along the chromosomes together when hmmlearn calls
        iter_from_X_lengths()
        Thus, X is combined binomial counts (col 1), size (col2), and
        probabilities (col 3 + 4).
        xs and ns have shape (n_obs in sequence, n_features==1).
        ps is the emission probabilities and has shape
            (n_components [states] in HMM == 2, n_features==1)."""
        assert type(X).__module__ == "numpy"
        xs = X[:, 0]
        ns = X[:, 1]
        ps = X[:, 2:4]
        ref = binom.logpmf(xs, ns, ps[:, 0])
        nonref = binom.logpmf(xs, ns, ps[:, 1])
        return np.matrix([ref, nonref]).T
