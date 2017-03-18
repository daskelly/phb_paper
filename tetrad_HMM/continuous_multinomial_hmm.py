"""Extend hmmlearn with an HMM that emits multinomially distributed
data that is conditional on total coverage at the site."""
import numpy as np
import sys
from hmmlearn import base
from hmmlearn.utils import log_normalize, normalize, iter_from_X_lengths
from sklearn.utils import check_array, check_random_state
from scipy.special import gammaln

# X is the data. Dim ~100k*4. rows = sites, cols = coverage
def log_dmultinomial(xs, ps):
    """Return the log of the multinomial probability density.
    xs has shape (n_obs in sequence, n_features).
    ps is the emission probabilities and has shape
        (n_components [states] in HMM, n_features)."""
    assert type(xs).__module__ == "numpy"
    assert type(ps).__module__ == "numpy"
    assert xs.shape[1] == ps.shape[1]
    # xs.shape[0] is possibly very big. Make sure we vectorize over it.
    n = xs.sum(axis=1)  # sum over rows
                        # i.e. sample size for each observation in sequence
    part1 = gammaln(n + 1)                                      # shape (n_obs,)
    part2 = np.apply_along_axis(gammaln, 0, xs + 1).sum(axis=1) # shape (n_obs,)
    part3 = np.dot(xs, np.log(ps.T)).T             # shape (n_components, n_obs)
    # for any given state, part3 is (xs * np.log(probs)).sum(axis=1)
    return (part1 - part2 + part3).T

# assert all(emissionp.sum(axis=1) == 1)
class ContinuousMultinomialHMM(base._BaseHMM):
    """Extends base._BaseHMM in the hmmlearn package for use with data
    where the emissions are multinomial data conditional on a given
    sample size at each site (sequence coverage).
    Significant basing of the code below on hmm.MultinomialHMM.

    est_transition_indices is for when you want to have a transition
    matrix that is partially estimated and partially fixed."""
    def __init__(self, n_components=1,
                 startprob_prior=1.0, transmat_prior=1.0,
                 algorithm="viterbi", random_state=None,
                 n_iter=10, tol=1e-2, verbose=False,
                 params="ste", init_params="ste",
                 simulated_sample_coverage=100,
                 est_transition_indices=None):
        base._BaseHMM.__init__(self, n_components,
                          startprob_prior=startprob_prior,
                          transmat_prior=transmat_prior,
                          algorithm=algorithm,
                          random_state=random_state,
                          n_iter=n_iter, tol=tol, verbose=verbose,
                          params=params, init_params=init_params)
        self.simulated_sample_coverage = simulated_sample_coverage
        if est_transition_indices is None:
            est_transition_indices = np.arange(n_components)
        self.est_transition_indices = est_transition_indices
        self.fix_transition_indices = np.array(list( \
            set(np.arange(n_components)) - set(est_transition_indices)))

    def _accumulate_sufficient_statistics(self, stats, X, framelogprob,
                                          posteriors, fwdlattice, bwdlattice):
        """Code called in hmmlearn.base._BaseHMM.fit()
        """
        old_t = np.copy(self.transmat_)
        super(ContinuousMultinomialHMM, self)._accumulate_sufficient_statistics(
            stats, X, framelogprob, posteriors, fwdlattice, bwdlattice)
        f = self.fix_transition_indices
        if len(f) > 0:
            # I'm running into this problem where I have an LOH state and
            # I want to fix transition probabilities so it is a sink state
            # But the function above sets them all to zero. So below is a hack
            # to get them to what they should be ...
            self.transmat_[f[:, None], f] = old_t[f[:, None], f]
        if 'e' in self.params:
            # stats['obs'] is a matrix with shape (n_components, n_features)
            # can think of in comparison to matrix of emission probabilities
            # which is about to be updated
            #
            # The dot product of posteriors and X will give
            # the un-normalized updated emission probability matrix.
            stats['obs'] += np.dot(posteriors.T, X)

    def _check(self):
        super(ContinuousMultinomialHMM, self)._check()

        self.emissionprob_ = np.atleast_2d(self.emissionprob_)
        n_features = getattr(self, "n_features", self.emissionprob_.shape[1])
        if self.emissionprob_.shape != (self.n_components, n_features):
            raise ValueError(
                "emissionprob_ must have shape (n_components, n_features)")
        else:
            self.n_features = n_features

    def _compute_log_likelihood(self, X):
        """Needs to return a matrix that is of shape
        (n_obs in sequence, n_components).
        For each component (state) it returns the log density. It does
        this for each observation (or set of observations of length n_features)
        in the sequence."""
        return log_dmultinomial(X, self.emissionprob_)

    def _do_mstep(self, stats):
        "customized for combination fixed and estimated transition probs"
        assert not 's' in self.params
        # update transition probs if specified in params
        # follows base._do_mstep()
        if 't' in self.params:
            t_ind = self.est_transition_indices
            cp = np.copy(self.transmat_[t_ind[:, None], t_ind])
            cp = np.where(cp == 0.0, cp, stats['trans'][t_ind[:, None], t_ind])
            normalize(cp, axis=1)
            self.transmat_[t_ind[:, None], t_ind] = cp
        # update emission probs if specified in params
        if 'e' in self.params:
            # the below code just normalizes the stats['obs'] matrix
            # so the rows sum to 1. The rows themselves represent the
            # probability of drawing from each of the n_features
            # features, for each component (of which there are n_components)
            sums = stats['obs'].sum(1)
            stats['obs'][sums == 0, :] = [0.25, 0.25, 0.25, 0.25]
            self.emissionprob_ = (stats['obs']
                                  / stats['obs'].sum(1)[:, np.newaxis])
        if np.any(np.isnan(self.transmat_)) or \
            np.any(np.isnan(self.emissionprob_)):
            print >> sys.stderr, stats['obs']
            raise Exception("Found nans")


    def _generate_sample_from_state(self, state, random_state=None):
        """Generate a random sample from a state. This is used in the
        hmmlearn.base._BaseHMM.sample() method."""
        N = self.simulated_sample_coverage
        eprobs = self.emissionprob_[state, :]
        return np.random.multinomial(N, eprobs, size=1)[0]

    def _init(self, X, lengths=None):
        """estimate emission probs.
        FYI: lengths is for when you have multiple sequences. E.g.
        sequences for different chroms that need to be run through the HMM
        separately."""
        super(ContinuousMultinomialHMM, self)._init(X, lengths=lengths)
        self.random_state = check_random_state(self.random_state)

        if 'e' in self.init_params:
            if not hasattr(self, "n_features"):
                self.n_features = X.shape[1]

            # code from hmm.MultinomialHMM._init
            # random initialization of emission probs
            self.emissionprob_ = self.random_state \
                .rand(self.n_components, self.n_features)
            # from each component, the emission probs should sum to 1:
            normalize(self.emissionprob_, axis=1)

    def _initialize_sufficient_statistics(self):
        """Sufficient stat for multinomial is the observations.
        Code cannibalized from hmm.MultinomialHMM"""
        stats = super(ContinuousMultinomialHMM, \
            self)._initialize_sufficient_statistics()
        stats['obs'] = np.zeros((self.n_components, self.n_features))
        return stats
