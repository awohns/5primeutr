# MIT License
#
# Copyright (c) 2024 Wilder Wohns
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
Infer the relationship between function and fitness of genetic variants
"""
import collections
import numpy as np
import pickle
import scipy.stats as stats

from tqdm import tqdm


class Likelihoods:
    """
    A class to manage the likelihoods of variant frequency and selection
    coefficients. These likelihoods are precomputed using the discrete-time
    Wright-Fisher model (DTWF) from Spence et al. (2023) as implemented in
    https://github.com/jeffspence/fastDTWF.
    This class provides an interface to load, subsample, normalize, and access
    the likelihood data computed by DTWF.
    """

    def __init__(self, likelihoods_path, S_GRID):
        self.likelihoods_path = likelihoods_path
        self.S_GRID = S_GRID

    def load(self, normalize=True):
        """
        Loads the likelihoods pickle file computed using dtwf.
        """
        if type(self.likelihoods_path) is dict:
            self.likelihoods = {}
            for sample_size, path in self.likelihoods_path.items():
                with open(path, "rb") as handle:
                    self.likelihoods[sample_size] = pickle.load(handle)
                if normalize:
                    lls_dict_normalized = {}
                    for key, val in tqdm(self.likelihoods[sample_size].items()):
                        # for col in range(val.shape[1]):
                        # Remove the row of frequency = 0
                        ll_curves_modified = val[1:-1, :-1]

                        # Calculate the sum of each column
                        column_sums = ll_curves_modified.sum(axis=0)

                        # Normalize each column so they sum to 1
                        ll_curves_modified = ll_curves_modified / column_sums
                        lls_dict_normalized[key] = ll_curves_modified
                    self.likelihoods[sample_size] = lls_dict_normalized
            # Check the shape of the S_GRID is the same as the columns of the likelihoods
            for sample_size, ll_dict in self.likelihoods.items():
                for key, ll_val in ll_dict.items():
                    assert ll_val.shape[1] == self.S_GRID.shape[0]
        else:
            with open(self.likelihoods_path, "rb") as handle:
                self.likelihoods = pickle.load(handle)
            if normalize:
                lls_dict_normalized = {}
                for key, val in tqdm(self.likelihoods.items()):
                    # for col in range(val.shape[1]):
                    # Remove the row of frequency = 0
                    ll_curves_modified = val[1:-1, :-1]

                    # Calculate the sum of each column
                    column_sums = ll_curves_modified.sum(axis=0)

                    # Normalize each column so they sum to 1
                    ll_curves_modified = ll_curves_modified / column_sums
                    lls_dict_normalized[key] = ll_curves_modified
                self.likelihoods = lls_dict_normalized

            # Check the shape of the S_GRID is the same as the columns of the likelihoods
            for key, ll_val in self.likelihoods.items():
                assert ll_val.shape[1] == self.S_GRID.shape[0]

    def compute_ll(
        self,
        c,
        sample_sizes,
        cur_betas,
        cur_sds,
        cur_s_het,
        cur_q,
        contexts,
        neutral_prob=0,
        progress=False,
        normalized=False,
    ):
        assert np.all(
            np.isclose(
                cur_betas.shape[0],
                np.array(
                    [
                        cur_betas.shape[0],
                        cur_sds.shape[0],
                        cur_s_het.shape[0],
                        cur_q.shape[0],
                        contexts.shape[0],
                    ]
                ),
            )
        )

        # make a list of the sample sizes in the likelihoods
        sample_size_keys = np.array(list(self.likelihoods.keys()))

        # We need to renormalize the likelihoods by dividing by all nonzero columns
        cache_denominators = {}
        for context in np.unique(contexts):
            for sample_size in sample_size_keys:
                cache_denominators[context + "_" + str(sample_size)] = np.zeros_like(
                    self.S_GRID
                )
                for closest_s, s_cur in enumerate(self.S_GRID):
                    cache_denominators[context + "_" + str(sample_size)][
                        closest_s
                    ] = np.sum(self.likelihoods[sample_size][context][1:, closest_s])

        cache_beta_range = {}
        for s_het in np.unique(cur_s_het):
            # We want to find values of beta that are equivalent to each value in the S_GRID
            # Given c and shet
            beta_range = np.sqrt(self.S_GRID / (s_het * c))
            beta_range = np.concatenate([beta_range[::2], -beta_range[1::2]])
            beta_range = beta_range[
                np.argsort(np.abs(beta_range))
            ]  # beta_range[np.argsort(np.abs(beta_range))]**2*1*cur_s_het[0]

            cache_beta_range[s_het] = beta_range

        total_ll = 0

        for beta, sd, s_het, cur_q, context, sample_sizes in tqdm(
            zip(cur_betas, cur_sds, cur_s_het, cur_q, contexts, sample_sizes),
            disable=not progress,
        ):
            if normalized:
                reference_allele_frequences = np.arange(1, sample_size + 1)
            else:
                reference_allele_frequences = np.arange(0, sample_size + 1)

            closest_sample_size = sample_size_keys[
                np.argmin(np.abs(sample_size_keys - sample_size))
            ]
            beta_range = cache_beta_range[s_het]
            # Calculate weights for each value of beta
            if sd != 0:
                weights = stats.norm(beta, sd).pdf(beta_range)
                if np.all(weights == 0):
                    # if all weights are 0 and the beta is greater than the beta_range
                    # assign all weight to highest beta_range value
                    if beta > np.max(beta_range):
                        weights[np.argmax(beta_range)] = 1
                    elif beta < np.min(beta_range):
                        weights[np.argmin(beta_range)] = 1
                    else:
                        print(
                            "Failing beta",
                            beta,
                            np.max(beta_range),
                            weights,
                            beta_range,
                        )
                        weights[np.argmin(np.abs(beta_range - beta))] = 1
                        if np.all(weights == 0):
                            print(
                                "Failing beta",
                                beta,
                                np.max(beta_range),
                                weights,
                                beta_range,
                            )
                            raise ValueError

                # Normalize weights to sum to 1
                weights = weights / np.sum(weights)
                assert ~np.any(np.isnan(weights)), weights
            else:
                weights = np.zeros_like(beta_range)
                weights[np.argmin(np.abs(beta - beta_range))] = 1

            # Find the index of the closest frequency in the reference list
            q_index = np.where(cur_q == reference_allele_frequences)[0][0]

            prob = 0
            denominator = 0

            for closest_s, s_cur in enumerate(self.S_GRID):
                if (
                    weights[closest_s] != 0
                    and self.likelihoods[closest_sample_size][context][
                        q_index, closest_s
                    ]
                    != 0
                ):
                    prob += np.exp(
                        np.log(
                            self.likelihoods[closest_sample_size][context][
                                q_index, closest_s
                            ]
                        )
                        + np.log(weights[closest_s])
                    )
                    denominator += np.exp(
                        np.log(
                            cache_denominators[context + "_" + str(sample_size)][
                                closest_s
                            ]
                        )
                        + np.log(weights[closest_s])
                    )

            if neutral_prob != 0:
                neutral_ll = self.likelihoods[closest_sample_size][context][q_index, 0]
                # Normalization: we want the unnormalized probability divided by the normalization
                numerator = (prob * (1 - neutral_prob)) + (neutral_ll * neutral_prob)
                denominator = (denominator * (1 - neutral_prob)) + (
                    cache_denominators[context + "_" + str(sample_size)][0]
                    * neutral_prob
                )
                # The last column (highest s value) has all of its probability mass on singletons.
                # This means we can have an issue where the probability is 0 and denominator is 0 if the frequency is non-singleton
                # Hacky fix here is to just add a very small probability
                if denominator == 0 and prob == 0:
                    total_ll += -np.log(np.nextafter(0, 1))
                else:
                    total_ll += -np.log(numerator / denominator + 1e-8)
            else:
                if normalized:
                    total_ll += -np.log(prob + 1e-8)
                else:
                    # The last column (highest s value) has all of its probability mass on singletons.
                    # This means we can have an issue where the probability is 0 and denominator is 0 if the frequency is non-singleton
                    # Hacky fix here is to just add a very small probability
                    if denominator == 0 and prob == 0:
                        total_ll += -np.log(np.nextafter(0, 1))
                    else:
                        total_ll += -np.log(prob / denominator + 1e-8)

        return total_ll

    def maximum_likelihood(
        self,
        betas,
        s_het,
        qs,
        contexts,
        sample_sizes,
        sds=None,
        neutral_prob=0,
        num_cs=15,
        max_cs=5,
        progress=False,
        normalized=False,
    ):
        cs = np.linspace(0.01, max_cs, num_cs)

        if type(neutral_prob) is np.ndarray:
            lls = collections.defaultdict(list)
            for n_prob in tqdm(
                neutral_prob, disable=not progress, desc="Neutral Probabilities"
            ):
                for c in tqdm(cs, disable=not progress, desc="Values of c"):
                    ll = self.compute_ll(
                        c,
                        sample_sizes,
                        betas,
                        sds,
                        s_het,
                        qs,
                        contexts,
                        neutral_prob=n_prob,
                        progress=False,
                        normalized=normalized,
                    )
                    lls[n_prob].append(ll)
        elif np.any([type(neutral_prob) is float, neutral_prob == 0]):
            lls = []
            for c in tqdm(cs, disable=not progress):
                ll = self.compute_ll(
                    c,
                    sample_sizes,
                    betas,
                    sds,
                    s_het,
                    qs,
                    contexts,
                    neutral_prob=neutral_prob,
                    progress=False,
                    normalized=normalized,
                )
                lls.append(ll)
        else:
            raise ValueError(
                "neutral_prob must be 0, a float or a numpy array of neutral probabilities to assess"
            )
        return cs, lls
