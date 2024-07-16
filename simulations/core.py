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
import pandas as pd
import pickle
import re
import scipy.stats as stats
import math
import warnings

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

    def __init__(self, likelihoods_path, S_GRID, mutation_rates_path):
        self.likelihoods_path = likelihoods_path
        self.S_GRID = S_GRID
        self.mutation_rates_df = pd.read_csv(mutation_rates_path, delimiter="\t")
        self.mutation_rates_df.loc[:, "full_context_nochrom"] = (
            self.mutation_rates_df["context"].str[0]
            + "("
            + self.mutation_rates_df["ref"]
            + ">"
            + self.mutation_rates_df["alt"]
            + ")"
            + self.mutation_rates_df["context"].str[-1]
            + "_"
            + self.mutation_rates_df["methylation_level"].astype(str)
        )

    # Function to parse the full context
    def parse_context(self, full_context):
        match = re.match(
            r"([ATGC])\(([ATGC])>([ATGC])\)([ATGC])_(\d+)([XA])", full_context
        )
        if match:
            context = match.group(1) + match.group(2) + match.group(4)
            ref = match.group(2)
            alt = match.group(3)
            methylation_level = int(match.group(5))
            return context, ref, alt, methylation_level, match.group(6), full_context
        else:
            raise ValueError("Context does not match expected pattern")

    def check_consistency(self):
        """
        Determine if likelihood dictionaries have the expected shape.
        Extract mutation rates for the full contexts.
        Create separate DataFrames for autosomal and X chromosomal contexts.
        """

        def extract_and_check_contexts(ll_dict, first_full_contexts, second_dim_shape):
            for key, ll_val in ll_dict.items():
                # Check the second dimension of the shape
                if second_dim_shape is None:
                    second_dim_shape = ll_val.shape[1]
                else:
                    assert (
                        ll_val.shape[1] == second_dim_shape
                    ), f"Inconsistent second dimension found: {ll_val.shape[1]} != {second_dim_shape}"

                # Check the shape of the S_GRID is the same as the second dimension of the likelihoods
                assert (
                    ll_val.shape[1] == self.S_GRID.shape[0]
                ), "Shape mismatch with S_GRID"

                # Extract full contexts
                ll_contexts = np.array(list(ll_dict.keys()))

                # If this is the first iteration, store the contexts
                if first_full_contexts is None:
                    first_full_contexts = ll_contexts
                else:
                    # Check consistency with the first encountered contexts
                    assert np.array_equal(
                        ll_contexts, first_full_contexts
                    ), "Inconsistent full contexts found"

            return first_full_contexts, second_dim_shape

        first_full_contexts = None
        second_dim_shape = None

        if all(isinstance(val, dict) for val in self.likelihoods.values()):
            # Handle the nested dictionary case
            for sample_size, ll_dict in self.likelihoods.items():
                first_full_contexts, second_dim_shape = extract_and_check_contexts(
                    ll_dict, first_full_contexts, second_dim_shape
                )
        else:
            # Handle the single dictionary case
            first_full_contexts, second_dim_shape = extract_and_check_contexts(
                self.likelihoods, first_full_contexts, second_dim_shape
            )

        # After the loop, set the final contexts and mutation_rates_ll
        self.ll_contexts = np.array(first_full_contexts)

        # Parse the contexts
        parsed_contexts = [self.parse_context(fc) for fc in first_full_contexts]

        # Filter the mutation rates DataFrame
        ll_mutation_rates_df = pd.DataFrame()

        for (
            context,
            ref,
            alt,
            methylation_level,
            chrom,
            full_context,
        ) in parsed_contexts:
            match_df = self.mutation_rates_df[
                (self.mutation_rates_df.context == context)
                & (self.mutation_rates_df.ref == ref)
                & (self.mutation_rates_df.alt == alt)
                & (self.mutation_rates_df.methylation_level == methylation_level)
            ]
            match_df = match_df.copy()  # Avoid SettingWithCopyWarning
            match_df.loc[:, "full_context"] = full_context

            ll_mutation_rates_df = pd.concat([ll_mutation_rates_df, match_df])
        contains_X = ll_mutation_rates_df["full_context"].str.contains("X")
        self.ll_mutation_rates_df_autosomes = ll_mutation_rates_df[
            ~contains_X
        ].sort_values(by="mu_snp")
        self.ll_mutation_rates_df_X = ll_mutation_rates_df[contains_X].sort_values(
            by="mu_snp"
        )
        # make a list of the sample sizes in the likelihoods
        self.sample_size_keys = np.array(list(self.likelihoods.keys()))

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
            self.check_consistency()
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

            self.check_consistency()

    def interpolate_ll(self, context, mu_snp, closest_sample_size):
        """
        If the trinucleotide context of a variant does not appear in the likelihood table,
        choose the closest two contexts and interpolate.
        Mutation rates differ for the autosomes and X chromosome. Only interpolate within the appropriate
        category for each variant.
        """
        # If the context for the variant is in the likelihoods there's no need for interpolation
        if context in self.ll_contexts:
            return self.likelihoods[closest_sample_size][context]

        if "X" in context:
            mutation_df = self.ll_mutation_rates_df_X
        else:
            mutation_df = self.ll_mutation_rates_df_autosomes
        if mu_snp < np.min(mutation_df["mu_snp"]):
            closest_context = mutation_df.loc[
                mutation_df["mu_snp"].idxmax(), "full_context"
            ]
            return self.likelihoods[closest_sample_size][closest_context]
        if mu_snp > np.max(mutation_df["mu_snp"]):
            closest_context = mutation_df.loc[
                mutation_df["mu_snp"].idxmax(), "full_context"
            ]
            return self.likelihoods[closest_sample_size][closest_context]

        # Compute distances and find the two closest contexts
        distances = np.abs(mu_snp - mutation_df["mu_snp"].values)
        sorted_indices = np.argsort(distances)
        context1 = mutation_df.iloc[sorted_indices[0]]["full_context"]
        context2 = mutation_df.iloc[sorted_indices[1]]["full_context"]

        # Compute weights for interpolation
        dist1 = distances[sorted_indices[0]]
        dist2 = distances[sorted_indices[1]]
        weight1 = 1 / dist1
        weight2 = 1 / dist2
        total_weight = weight1 + weight2

        # Interpolate the likelihood curves
        curve1 = self.likelihoods[closest_sample_size][context1]
        curve2 = self.likelihoods[closest_sample_size][context2]
        interpolated_curve = (curve1 * weight1 + curve2 * weight2) / total_weight

        return interpolated_curve

    def closest_ll(self, context, mu_snp, closest_sample_size):
        """
        If the trinucleotide context of a variant does not appear in the likelihood table,
        choose the closest contexts.
        Mutation rates differ for the autosomes and X chromosome. Only choose the appropriate
        category for each variant.
        """
        # If the context for the variant is in the likelihoods there's no need for interpolation
        if context in self.ll_contexts:
            return self.likelihoods[closest_sample_size][context]

        if "X" in context:
            mutation_df = self.ll_mutation_rates_df_X
        else:
            mutation_df = self.ll_mutation_rates_df_autosomes
        if mu_snp < np.min(mutation_df["mu_snp"]):
            closest_context = mutation_df.loc[
                mutation_df["mu_snp"].idxmax(), "full_context"
            ]
            return self.likelihoods[closest_sample_size][closest_context]
        if mu_snp > np.max(mutation_df["mu_snp"]):
            closest_context = mutation_df.loc[
                mutation_df["mu_snp"].idxmax(), "full_context"
            ]
            return self.likelihoods[closest_sample_size][closest_context]

        # Compute distances and find the two closest contexts
        distances = np.abs(mu_snp - mutation_df["mu_snp"].values)
        sorted_indices = np.argsort(distances)
        context = mutation_df.iloc[sorted_indices[0]]["full_context"]
        curve = self.likelihoods[closest_sample_size][context]

        return curve

    def create_ll_cache(self, contexts, cur_q):
        # We need to renormalize the likelihoods by dividing by all nonzero columns
        self.cache_denominators = {}
        for context in tqdm(np.unique(contexts)):
            for sample_size in self.sample_size_keys:
                if sample_size >= np.max(cur_q):
                    closest_sample_size = self.sample_size_keys[-1]
                else:
                    closest_sample_size = self.sample_size_keys[
                        np.searchsorted(
                            self.sample_size_keys,
                            [
                                sample_size,
                            ],
                            side="right",
                        )[0]
                    ]

                self.cache_denominators[
                    context + "_" + str(closest_sample_size)
                ] = np.zeros_like(self.S_GRID)
                mut_rate = self.mutation_rates_df[
                    self.mutation_rates_df.full_context_nochrom == context[:-1]
                ]["mu_snp"]
                assert mut_rate.shape[0] == 1
                closest_ll = self.closest_ll(
                    context, mut_rate.iloc[0], closest_sample_size
                )
                for closest_s, s_cur in enumerate(self.S_GRID):
                    self.cache_denominators[context + "_" + str(closest_sample_size)][
                        closest_s
                    ] = np.sum(closest_ll[1:, closest_s])

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
        self.create_ll_cache(contexts, cur_q)

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

        for beta, sd, s_het, cur_q, context, sample_size in tqdm(
            zip(cur_betas, cur_sds, cur_s_het, cur_q, contexts, sample_sizes),
            disable=not progress,
        ):
            if normalized:
                reference_allele_frequences = np.arange(1, sample_size + 1)
            else:
                reference_allele_frequences = np.arange(0, sample_size + 1)

            if sample_size >= np.max(cur_q):
                closest_sample_size = self.sample_size_keys[-1]
            else:
                closest_sample_size = self.sample_size_keys[
                    np.searchsorted(
                        self.sample_size_keys,
                        [
                            sample_size,
                        ],
                        side="right",
                    )[0]
                ]
            # print("closest", closest_sample_size, sample_size)
            beta_range = cache_beta_range[s_het]
            # Calculate weights for each value of beta
            if ~np.isclose(sd, 0):
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
                            sd,
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

                # avoid floating point errors
                if np.sum(weights) > 1e-300:
                    weights = weights / np.sum(weights)
                else:
                    # Check that all values are 0 except 1
                    assert np.sum(weights == 0) == (weights.shape[0] - 1), weights
                    assert weights[np.argmin(np.abs(beta_range - beta))] != 0
                    weights[np.argmin(np.abs(beta_range - beta))] = 1
                    print("avoided overflow")

                assert ~np.any(np.isnan(weights)), weights
            else:
                weights = np.zeros_like(beta_range)
                weights[np.argmin(np.abs(beta - beta_range))] = 1

            # Find the index of the closest frequency in the reference list
            q_index = np.where(cur_q == reference_allele_frequences)[0][0]

            prob = 0
            denominator = 0

            mut_rate = self.mutation_rates_df[
                self.mutation_rates_df.full_context_nochrom == context[:-1]
            ]["mu_snp"]
            assert mut_rate.shape[0] == 1
            cur_ll_curve = self.closest_ll(
                context, mut_rate.iloc[0], closest_sample_size
            )

            for closest_s, s_cur in enumerate(self.S_GRID):
                try:
                    cur_ll_curve[q_index, closest_s]
                except:
                    print(closest_sample_size, q_index, closest_s)
                if weights[closest_s] != 0 and cur_ll_curve[q_index, closest_s] != 0:
                    prob += np.exp(
                        np.log(cur_ll_curve[q_index, closest_s])
                        + np.log(weights[closest_s])
                    )
                    np.seterr(invalid="raise")
                    denominator += np.exp(
                        np.log(
                            self.cache_denominators[
                                context + "_" + str(closest_sample_size)
                            ][closest_s]
                        )
                        + np.log(weights[closest_s])
                    )

            if neutral_prob != 0:
                neutral_ll = cur_ll_curve[q_index, 0]
                # Normalization: we want the unnormalized probability divided by the normalization
                numerator = (prob * (1 - neutral_prob)) + (neutral_ll * neutral_prob)
                denominator = (denominator * (1 - neutral_prob)) + (
                    self.cache_denominators[context + "_" + str(closest_sample_size)][0]
                    * neutral_prob
                )
                # The last column (highest s value) has all of its probability mass on singletons.
                # This means we can have an issue where the probability is 0 and denominator is 0 if the frequency is non-singleton
                # Hacky fix here is to just add a very small probability
                if denominator == 0 and prob == 0:
                    raise ValueError
                    # total_ll += -np.log(np.nextafter(0, 1))
                else:
                    total_ll += -np.log(numerator / denominator + 1e-8)
                    print(-np.log(numerator / denominator + 1e-8))
            else:
                if normalized:
                    total_ll += -np.log(prob + 1e-8)
                else:
                    # The last column (highest s value) has all of its probability mass on singletons.
                    # This means we can have an issue where the probability is 0 and denominator is 0 if the frequency is non-singleton
                    # Hacky fix here is to just add a very small probability
                    if denominator == 0 and prob == 0:
                        continue

                    else:
                        total_ll += -np.log(prob / denominator + 1e-8)

        return total_ll

    def likelihood_ratio_test(self, neutral_ll, alt_ll):
        # Likelihood ratio test is -2 * ln(L_0) - ln(l_alt)
        chi2_val = 2 * (neutral_ll - alt_ll)
        return stats.chi2.sf(chi2_val, 1)

    def fisher_information(self, lls, mle, cs):
        print(lls, mle[0])
        diffs = np.diff(lls[mle[0]])
        min_ll_index = np.argmin(lls[mle[0]][mle[1]])
        return np.sqrt(
            1
            / (
                np.abs(diffs[min_ll_index] - diffs[min_ll_index - 1])
                / (np.diff(cs)[0] ** 2)
            )
        )

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
        compute_fisher_information=False,
    ):
        cs = np.linspace(0.01, max_cs, num_cs)
        # cs = np.array([1.6])

        if type(neutral_prob) is np.ndarray:
            lls = collections.defaultdict(list)
            mle = (neutral_prob[0], cs[0], 1e100)
            for n_prob in tqdm(
                neutral_prob, disable=not progress, desc="Neutral Probabilities"
            ):
                for c_index, c in enumerate(cs):
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
                    if ll < mle[2]:
                        mle = (n_prob, c_index, ll)

            print(
                "Likelihood ratio test. Alternative: n_prob: {}, c: {}, p-value {}".format(
                    mle[0],
                    cs[mle[1]],
                    self.likelihood_ratio_test(lls[0][0], lls[mle[0]][mle[1]]),
                )
            )
            print(
                "Likelihood ratio test. Alternative: n_prob: {}, c: {}, p-value {}".format(
                    mle[0],
                    cs[mle[1]],
                    self.likelihood_ratio_test(
                        lls[mle[0]][np.argmin(np.abs(cs - 1))], lls[mle[0]][mle[1]]
                    ),
                )
            )

            ll_ratio_neutral = (
                mle[0],
                cs[mle[1]],
                self.likelihood_ratio_test(lls[0][0], lls[mle[0]][mle[1]]),
            )
            ll_ratio_c_1 = self.likelihood_ratio_test(
                lls[mle[0]][np.argmin(np.abs(cs - 1))], lls[mle[0]][mle[1]]
            )
            if compute_fisher_information:
                print(
                    "MLE Estimate c: {}, Margin of error: {}".format(
                        cs[mle[1]], self.fisher_information(lls, mle, cs)
                    )
                )
                fisher_information = self.fisher_information(lls, mle, cs)
            else:
                fisher_information = None
            return cs, lls, mle, ll_ratio_neutral, ll_ratio_c_1, fisher_information
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
            min_c = np.argmin(lls)
            fisher_information = self.fisher_information(
                [lls], (0, min_c, np.min(cs)), cs
            )
            return cs, lls, fisher_information
        else:
            raise ValueError(
                "neutral_prob must be 0, a float or a numpy array of neutral probabilities to assess"
            )
