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
from scipy.stats import norm

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

    def load(self, progress=False, normalize=True):
        """
        Loads the likelihoods pickle file computed using dtwf.
        """
        self.normalized = normalize
        if type(self.likelihoods_path) is dict:
            self.likelihoods = {}
            self.normalization_constants = {}
            for sample_size, path in self.likelihoods_path.items():
                with open(path, "rb") as handle:
                    self.likelihoods[sample_size] = pickle.load(handle)
                if normalize:
                    lls_dict_normalized = {}
                    for key, val in tqdm(
                        self.likelihoods[sample_size].items(), disable=not progress
                    ):
                        # for col in range(val.shape[1]):
                        # Remove the row of frequency = 0
                        ll_curves_modified = val[1:-1, :-1]

                        # Calculate the sum of each column
                        column_sums = ll_curves_modified.sum(axis=0)

                        # Normalize each column so they sum to 1
                        ll_curves_modified = ll_curves_modified / column_sums
                        lls_dict_normalized[key] = ll_curves_modified
                    self.likelihoods[sample_size] = lls_dict_normalized
                else:
                    normalization_constants = {}
                    for key, val in tqdm(
                        self.likelihoods[sample_size].items(), disable=not progress
                    ):
                        # for col in range(val.shape[1]):
                        # Remove the row of frequency = 0
                        ll_curves_modified = val[1:, :]

                        # Calculate the sum of each column
                        column_sums = ll_curves_modified.sum(axis=0)

                        normalization_constants[key] = column_sums
                    self.normalization_constants[sample_size] = normalization_constants
                # else:
                #     # Remove last column
                #     lls_dict_remove_last_col = {}
                #     for key, val in tqdm(
                #         self.likelihoods[sample_size].items(), disable=not progress
                #     ):  
                #         ll_curves_modified = val[0:-1, :-1]

                #         # Calculate the sum of each column
                #         column_sums = ll_curves_modified.sum(axis=0)

                #         # Normalize each column so they sum to 1
                #         ll_curves_modified = ll_curves_modified / column_sums
                #         lls_dict_remove_last_col[key] = ll_curves_modified
                #     self.likelihoods[sample_size] = lls_dict_remove_last_col
            self.check_consistency()
        else:
            with open(self.likelihoods_path, "rb") as handle:
                self.likelihoods = pickle.load(handle)
            if normalize:
                lls_dict_normalized = {}
                for key, val in tqdm(self.likelihoods.items(), disable=not progress):
                    # for col in range(val.shape[1]):
                    # Remove the row of frequency = 0
                    ll_curves_modified = val[1:-1, :-1]

                    # Calculate the sum of each column
                    column_sums = ll_curves_modified.sum(axis=0)

                    # Normalize each column so they sum to 1
                    ll_curves_modified = ll_curves_modified / column_sums
                    lls_dict_normalized[key] = ll_curves_modified
                self.likelihoods = lls_dict_normalized
            else:
                normalization_constants = {}
                for key, val in tqdm(
                    self.likelihoods.items(), disable=not progress
                ):
                    # for col in range(val.shape[1]):
                    # Remove the row of frequency = 0
                    ll_curves_modified = val[1:, :]

                    # Calculate the sum of each column
                    column_sums = ll_curves_modified.sum(axis=0)

                    normalization_constants[key] = column_sums
                self.normalization_constants = normalization_constants
            # else:
            #     lls_dict_remove_last_col = {}
            #     for key, val in tqdm(self.likelihoods.items(), disable=not progress):
            #         ll_curves_modified = val[0:-1, :-1]

            #         # Calculate the sum of each column
            #         column_sums = ll_curves_modified.sum(axis=0)

            #         # Normalize each column so they sum to 1
            #         ll_curves_modified = ll_curves_modified / column_sums
            #         lls_dict_remove_last_col[key] = ll_curves_modified
            #     self.likelihoods = lls_dict_remove_last_col

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
            return (closest_sample_size, context)

        if "X" in context:
            mutation_df = self.ll_mutation_rates_df_X
        else:
            mutation_df = self.ll_mutation_rates_df_autosomes
        if mu_snp < np.min(mutation_df["mu_snp"]):
            closest_context = mutation_df.loc[
                mutation_df["mu_snp"].idxmax(), "full_context"
            ]
            return (closest_sample_size, closest_context)
        if mu_snp > np.max(mutation_df["mu_snp"]):
            closest_context = mutation_df.loc[
                mutation_df["mu_snp"].idxmax(), "full_context"
            ]
            return (closest_sample_size, closest_context)

        # Compute distances and find the two closest contexts
        distances = np.abs(mu_snp - mutation_df["mu_snp"].values)
        sorted_indices = np.argsort(distances)
        context = mutation_df.iloc[sorted_indices[0]]["full_context"]
        # curve = self.likelihoods[closest_sample_size][context]

        return (closest_sample_size, context)

    def closest_context(self, context):
        """
        If the trinucleotide context of a variant does not appear in the likelihood table,
        choose the closest contexts.
        """
        if context in self.ll_contexts:
            return context
        if "X" in context:
            mutation_df = self.ll_mutation_rates_df_X
        else:
            mutation_df = self.ll_mutation_rates_df_autosomes
        mu_snp = self.mutation_rates_df[
            self.mutation_rates_df.full_context_nochrom == context[:-1]
        ]["mu_snp"].values
        assert mu_snp.shape[0] == 1
        if mu_snp < np.min(mutation_df["mu_snp"]):
            closest_context = mutation_df.loc[
                mutation_df["mu_snp"].idxmax(), "full_context"
            ]
            return closest_context
        if mu_snp > np.max(mutation_df["mu_snp"]):
            closest_context = mutation_df.loc[
                mutation_df["mu_snp"].idxmax(), "full_context"
            ]
            return closest_context

        # Compute distances and find the two closest contexts
        distances = np.abs(mu_snp - mutation_df["mu_snp"].values)
        sorted_indices = np.argsort(distances)
        closest_context = mutation_df.iloc[sorted_indices[0]]["full_context"]

        return closest_context

    def prior_s(self, c_value):
        self.all_beta_ranges = []
        self.all_weights = []

        # Precompute beta ranges for all unique combinations
        for s_het, beta in zip(self.s_het, self.betas):
            pos_beta_range = np.sqrt(self.S_GRID / (s_het * c_value))
            neg_beta_range = -np.sqrt(self.S_GRID / (s_het * c_value))
            self.all_beta_ranges.append(np.concatenate([np.flip(neg_beta_range), pos_beta_range]))
            # self.all_beta_ranges.append(np.sort(beta_range))
        self.all_beta_ranges = np.array(self.all_beta_ranges)
        loc = np.expand_dims(self.betas, axis=-1)  # Shape becomes (80330, 1)
        scale = np.expand_dims(self.sds, axis=-1)  # Shape becomes (80330, 1)
        self.all_weights = norm.pdf(self.all_beta_ranges, loc=loc, scale=scale)

        # For rows where SD is close to 0, assign all weight to closest value
        sd_zero = np.isclose(self.sds, 0)
        self.all_weights[sd_zero] = np.zeros_like(self.all_weights[sd_zero])
        expanded_betas = self.betas[sd_zero][:, np.newaxis]
        min_indices = np.argmin(np.abs(expanded_betas - self.all_beta_ranges[sd_zero]), axis=1)
        self.all_weights[sd_zero, min_indices] = 1

        assert ~np.any(np.isnan(self.all_weights))

    def create_ll_cache(self, closest_contexts, closest_sample_sizes, progress=False):
        # Precalculate P(S|AF>0) = P(S) * P(AF>0|S) / Constant
        # We need to renormalize the likelihoods by dividing by all nonzero columns
        self.cache_denominators = []
        
        for index, (closest_context, closest_sample_size) in tqdm(
            enumerate(zip(closest_contexts, closest_sample_sizes)), disable= not progress):
            closest_ll = self.likelihoods[int(closest_sample_size)][closest_context]
            # for closest_s, s_cur in enumerate(self.S_GRID):
            llhood = self.normalization_constants[closest_sample_size][closest_context]#np.sum(closest_ll[1:, :], axis=0)
            llhood = np.concatenate([np.flip(llhood), llhood])
            p_s = self.all_weights[index]
            numerator = llhood * p_s 
            self.cache_denominators.append(numerator / np.sum(numerator))

    def compute_posterior_s(self, neutral_prob):
        self.posterior_s = []
        self.posterior_s_0 = []
        for index, (closest_context, closest_sample_size, beta) in enumerate(
            zip(self.closest_contexts, self.closest_sample_sizes, self.betas)):
            beta_range = self.all_beta_ranges[index]
            beta_index = np.argmin(np.abs(beta_range - beta))
            af_greater_0 = self.normalization_constants[closest_sample_size][closest_context]
            af_greater_0 = np.concatenate([np.flip(af_greater_0), af_greater_0])
            posterior_s = (1 - neutral_prob) * af_greater_0[beta_index]
            posterior_s_0 = neutral_prob * self.normalization_constants[closest_sample_size][closest_context][0]
            constant = posterior_s + posterior_s_0
            if constant != 0:
                self.posterior_s.append(posterior_s/constant)
                self.posterior_s_0.append(posterior_s_0/constant)
            else:
                self.posterior_s.append(0)
                self.posterior_s_0.append(0)



    def load_data(
        self,
        sample_sizes,
        betas,
        sds,
        s_het,
        qs,
        contexts,
        progress=False,
    ):
        """
        Load the data required to perform inference.
        Assign contexts and sample sizes to their closest values in Likelihoods.
        """
        assert np.all(
            np.isclose(
                len(betas),
                np.array(
                    [
                        len(betas),
                        len(sds),
                        len(s_het),
                        len(qs),
                        len(contexts),
                    ]
                ),
            )
        )

        self.closest_sample_sizes = []
        self.closest_contexts = []
        for index, (context, sample_size) in enumerate(zip(contexts, sample_sizes)):
            # if sample_size >= np.max(qs):
            # closest_sample_size = self.sample_size_keys[-1]
            if sample_size in self.sample_size_keys:
                closest_sample_size = sample_size
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
            self.closest_sample_sizes.append(closest_sample_size)
            self.closest_contexts.append(self.closest_context(context))
        self.closest_contexts = np.array(self.closest_contexts)
        # self.create_ll_cache(
        #     self.closest_contexts, self.closest_sample_sizes, progress=progress
        # )
        self.sample_sizes = np.array(sample_sizes)
        self.closest_sample_sizes = np.array(self.closest_sample_sizes)
        self.betas = np.array(betas)
        self.sds = np.array(sds)
        self.s_het = np.array(s_het)
        self.qs = np.array(qs)

        if self.normalized:
            self.qs = self.qs.astype(int) - 1
        else:
            # if not normalized 
            self.qs = self.qs.astype(int)
        self.contexts = np.array(contexts)

        self.ll_indices = []
        self.param_history = []
        self.likelihood_history = []


    def quadratic_s(self, c, shet, beta):
        """
        Return the selection coefficient s, under a quadratic model where s = c * shet * beta**2
        """
        return c * shet * beta ** 2

    def compute_ll_piecewise_bin_no_se_no_neutral(self, c, bin_left, bin_right, progress=True):
        """
        Compute likelihood of data in a particular bin given c, without considering standard errors
        from ash or a neutral proportion. Bin edges are inclusive on the left and exclusive on the
        right, i.e.: [bin_left, bin_right)
        """
        assert self.normalized
        # determine which data points to include based on given bin
        bool_include = np.logical_and(self.betas >= bin_left, self.betas < bin_right)
        # calculate s given c and shet
        s = self.s_het[bool_include] * c
        # find the closest s index in the likelihood table for each included data point
        closest_s_index = [np.argmin(np.abs(self.S_GRID - s_val)) for s_val in s]
        # only include frequency, sample_size, and context data if they fall in the bin
        q_include = self.qs[bool_include]
        sample_size_include = self.closest_sample_sizes[bool_include]
        context_include = self.closest_contexts[bool_include]
        # iterate over data points, multiplying likelihoods to find the total_ll
        total_ll = 0
        for s_index, q_index, sample_size, context in zip(closest_s_index, q_include, sample_size_include, context_include):
            llhood = self.likelihoods[sample_size][context][q_index, s_index]
            total_ll += np.log(llhood)
        return - total_ll

    def compute_ll_piecewise_bin_no_se(self, c, neutral_prop, bin_left, bin_right, progress=True):
        """
        Compute likelihood of data in a particular bin given c, without considering standard errors
        from ash. Bin edges are inclusive on the left and exclusive on the
        right, i.e.: [bin_left, bin_right)
        """
        assert self.normalized
        # determine which data points to include based on given bin
        bool_include = np.logical_and(self.betas >= bin_left, self.betas < bin_right)
        # calculate s given c and shet
        s = self.s_het[bool_include] * c
        # find the closest s index in the likelihood table for each included data point
        closest_s_index = [np.argmin(np.abs(self.S_GRID - s_val)) for s_val in s]
        # only include frequency, sample_size, and context data if they fall in the bin
        q_include = self.qs[bool_include]
        sample_size_include = self.closest_sample_sizes[bool_include]
        context_include = self.closest_contexts[bool_include]
        # iterate over data points, multiplying likelihoods to find the total_ll
        total_ll = 0
        for s_index, q_index, sample_size, context in zip(closest_s_index, q_include, sample_size_include, context_include):
            llhood = self.likelihoods[sample_size][context][q_index, s_index]
            llhood_neutral = self.likelihoods[sample_size][context][q_index, 0]
            total_ll += np.log((llhood * (1 - neutral_prop)) + (llhood_neutral * neutral_prop))
        return - total_ll

    def compute_ll_quadratic(self, params, progress=True):
        """
        Compute likelihood under a quadratic model
        """
        results = []
        if len(params) == 3:
            c_positive = params[0]
            c_negative = params[1]
            neutral_prob = params[2]
            separate_cs = True
        elif len(params) == 2:
            c = params[0]
            neutral_prob = params[1]
            separate_cs = False
        elif len(params) == 1:
            c = params[0]
            neutral_prob = 0
            separate_cs = False
        self.prior_s(c) 
        if not self.normalized:
            self.create_ll_cache(self.closest_contexts, self.closest_sample_sizes)
        # # cache_beta_range = {}
        # # cache_weights = {}
        # self.all_beta_ranges = []
        # self.all_weights = []
        # # combined_array = np.column_stack((self.s_het, self.betas))
        # # unique_combinations = np.unique(combined_array, axis=0)

        # # Precompute beta ranges for all unique combinations
        # # for s_het, beta in tqdm(unique_combinations):
        # for s_het, beta in zip(self.s_het, self.betas):
        #     if separate_cs:
        #         c_value = c_positive if beta > 0 else c_negative
        #     else:
        #         c_value = c
        #     beta_range = np.sqrt(self.S_GRID / (s_het * c_value))
        #     self.all_beta_ranges.append(beta_range)
        # self.all_beta_ranges = np.array(self.all_beta_ranges)
        # loc = np.expand_dims(np.abs(self.betas), axis=-1)  # Shape becomes (80330, 1)
        # scale = np.expand_dims(self.sds, axis=-1)  # Shape becomes (80330, 1)
        # self.all_weights = norm.pdf(self.all_beta_ranges, loc=loc, scale=scale)

        # # CHECK THIS PART
        # # if all weights are 0 and the beta is greater than the beta_range
        # # assign all weight to highest beta_range value
        # zero_rows = np.all(self.all_weights == 0, axis=1)

        # # Expand dimensions to enable broadcasting
        # expanded_betas = self.betas[zero_rows][:, np.newaxis]

        # # Compute the difference and find the index of the minimum value
        # min_indices = np.argmin(np.abs(expanded_betas - self.all_beta_ranges[zero_rows]), axis=1)

        # # Assign weight to the closest beta range value
        # self.all_weights[zero_rows, min_indices] = 1

        # # For rows where SD is close to 0, assign all weight to closest value
        # sd_zero = np.isclose(self.sds, 0)
        # self.all_weights[sd_zero] = np.zeros_like(self.all_weights[sd_zero])
        # expanded_betas = self.betas[sd_zero][:, np.newaxis]
        # min_indices = np.argmin(np.abs(expanded_betas - self.all_beta_ranges[sd_zero]), axis=1)
        # self.all_weights[sd_zero, min_indices] = 1

        # # Sum along axis 1 and reshape for broadcasting
        # weights_sum = np.sum(self.all_weights, axis=1)[:, np.newaxis]
        
        # # Normalize weights
        # self.all_weights = self.all_weights / weights_sum
        # assert ~np.any(np.isnan(self.all_weights))

        # combined_array = np.column_stack((self.s_het, self.betas, self.sds))
        # unique_combinations = np.unique(combined_array, axis=0)
        # for s_het, beta, sd in tqdm(unique_combinations):
        #     beta_range = cache_beta_range[(s_het, beta)]

        #     if not np.isclose(sd, 0):
        #         # weights = stats.norm(beta, sd).pdf(beta_range)
        #         # weights = np.ones_like(beta_range)
        #         weights = norm.pdf(beta_range, loc=beta, scale=sd)
        #         if np.all(weights == 0):
        #             if beta > np.max(beta_range):
        #                 weights[np.argmax(beta_range)] = 1
        #             elif beta < np.min(beta_range):
        #                 weights[np.argmin(beta_range)] = 1
        #             else:
        #                 weights[np.argmin(np.abs(beta_range - beta))] = 1
        #                 if np.all(weights == 0):
        #                     raise ValueError("Failed to assign weights.")

        #         # if np.sum(weights) > 1e-300:
        #         #     weights /= np.sum(weights)
        #         # else:
        #         #     weights = np.zeros_like(beta_range)
        #         #     weights[np.argmin(np.abs(beta - beta_range))] = 1

        #         # assert not np.any(np.isnan(weights)), f"Weights contain NaN values: {weights}"
        #     else:
        #         weights = np.zeros_like(beta_range)
        #         weights[np.argmin(np.abs(beta - beta_range))] = 1

        #     cache_weights[(s_het, beta, sd)] = weights

        total_ll = 0

        # Optimization opportunity: iterate through groups of context/sample sizes
        # Then can look up all the data for each combo at once
        for index, (beta, sd, s_het, q_index, context, sample_size) in tqdm(
            enumerate(
                zip(
                    self.betas,
                    self.sds,
                    self.s_het,
                    self.qs,
                    self.closest_contexts,
                    self.closest_sample_sizes,
                )
            ),
            disable=not progress,
        ):

            # if not np.isclose(sd, 0):
            
                # ll_dict.cache_denominators contains P(AF>0|s)
    
                
                # prob = np.sum(
                #     np.exp(
                #         np.log(self.likelihoods[sample_size][context][q_index])
                #         + np.log(weights)
                #     )#[weights != 0]
                # )
                # denominator = np.sum(self.cache_denominators[context + "_" + str(sample_size)] * weights)
                # denominator = denominator / np.sum(denominator)
                # denominator = np.sum(
                    # np.exp(
                        # np.log(
                            # self.cache_denominators[context + "_" + str(sample_size)]
                        # )
                        # + np.log(weights)
                    # )#[weights != 0]
                # )

            # else:

                # prob = 0
                # denominator = 0
                # cur_ll_curve = self.likelihoods[sample_size][context]
                # beta_range = self.all_beta_ranges[index]  # [(s_het, beta)]
                # # print(np.argmin(np.abs((c * s_het * beta ** 2) - self.S_GRID)))
                # beta_index = np.argmin(np.abs(np.abs(beta) - beta_range))

                # prob = cur_ll_curve[q_index, beta_index]
                # # print(prob, sample_size, context, q_index, beta_index)                
                # denominator = self.cache_denominators[context + "_" + str(sample_size)][
                #     beta_index
                # ]

                # for closest_s, s_cur in enumerate(self.S_GRID):
                #     if weights[closest_s] != 0 and cur_ll_curve[q_index, closest_s] != 0:
                #         prob += np.exp(
                #             np.log(cur_ll_curve[q_index, closest_s])
                #             + np.log(weights[closest_s])
                #         )
                #         np.seterr(invalid="raise")
                #         denominator += np.exp(
                #             np.log(
                #                 self.cache_denominators[
                #                     context + "_" + str(closest_sample_size)
                #                 ][closest_s]
                #             )
                #             + np.log(weights[closest_s])
                #         )
            if not self.normalized:
                weights = self.all_weights[index]  # cache_weights[(s_het, beta, sd)]
                llhood = self.likelihoods[sample_size][context][q_index, :]
                # Note the likelihoods are determined with a grid of S, but S = c * beta **2 * shet, and we want it in terms of beta, so we have to 
                # account for positive and negative beta
                numerator = np.concatenate([np.flip(llhood), llhood])
                numerator = numerator[1:] #* np.diff(self.all_beta_ranges[index]) * self.all_weights[index][1:]
                # The denominator is p(AF>0|s)
                # denominator = np.sum(self.likelihoods[sample_size][context][1:, :], axis=0)
                denominator = self.normalization_constants[sample_size][context]
                denominator = np.concatenate([np.flip(denominator), denominator])
                denominator = denominator[1:] #* np.diff(self.all_beta_ranges[index]) * self.all_weights[index][1:]

                llhood = (numerator/denominator) #* np.diff(self.all_beta_ranges[index])# * self.all_weights[index][1:]
                # print(numerator, denominator)
                llhood[np.isnan(llhood)] = 0
                prob = np.nansum(
                                (llhood) *
                                 self.cache_denominators[index][1:])
            if neutral_prob != 0:
                neutral_ll = self.likelihoods[sample_size][context][q_index, 0]
                if self.normalized:
                    numerator = (prob * (1 - neutral_prob)) + (neutral_ll * neutral_prob)
                    total_ll += numerator
                    results.append(numerator)
                else:
                    # Normalization: we want the unnormalized probability divided by the normalization    
                    # Probability of s = 0 given AF > 0
                    beta_range = self.all_beta_ranges[index]
                    norm_constant = self.normalization_constants[sample_size][context]
                    posterior = ((self.all_weights[index][np.argmin(np.abs(beta_range-0))] * norm_constant[0]) /
                        np.sum(norm_constant))
                    neutral = (neutral_ll / norm_constant[0]) #* posterior
                    # denominator = (denominator * (1 - neutral_prob)) + (
                        # self.cache_denominators[context + "_" + str(sample_size)][0]
                        # * neutral_prob
                    # )
                    # The last column (highest s value) has all of its probability mass on singletons.
                    # This means we can have an issue where the probability is 0 and denominator is 0 if the frequency is non-singleton
                    # Hacky fix here is to just add a very small probability
                    # if denominator == 0 and prob == 0:
                        # raise ValueError
                        # total_ll += -np.log(np.nextafter(0, 1))
                    # else:
                    try:
                        if -np.log((prob * (1 - neutral_prob)) + neutral) == -np.inf:
                            # print(prob, neutral, neutral_ll, norm_constant, beta, beta_index)
                            continue
                        # print(-np.log(prob * (1 - neutral_prob)) + (neutral * neutral_prob), prob, neutral, neutral_prob)
                        total_ll += -np.log((prob * (1 - neutral_prob)) + (neutral * neutral_prob))
                        # total_ll += numerator / denominator
                        results.append(numerator / denominator)
                    except FloatingPointError as e:
                        print(f"Error encountered: {e}")
                        print(f"Numerator: {numerator}")
                        print(f"Denominator: {denominator}")
                        print(c, neutral_prob)
                        raise  # Re-raise the error after logging the information
            else:
                if self.normalized:
                    # beta_range = self.all_beta_ranges[index]  # [(s_het, beta)]
                    beta_range = np.sqrt(self.S_GRID / (s_het * c))
                    beta_index = np.argmin(np.abs(np.abs(beta) - beta_range))
                    llhood = self.likelihoods[sample_size][context][q_index, beta_index]
                    # Note the likelihoods are determined with a grid of S, but S = c * beta **2 * shet, and we want it in terms of beta, so we have to 
                    # account for positive and negative beta
                    # llhood = np.concatenate([np.flip(llhood), llhood])
                    # print(c, beta, beta_index, q_index, llhood)
                    total_ll += -np.log(llhood)  # + 1e-8)
                    # raise ValueError
                    # results.append(prob)
                else:
                    
                    # The last column (highest s value) has all of its probability mass on singletons.
                    # This means we can have an issue where the probability is 0 and denominator is 0 if the frequency is non-singleton
                    # Hacky fix here is to just add a very small probability
                    # if denominator == 0 and prob == 0:

                        # print(sample_size, context, beta_index, q_index, beta, sd, s_het, -np.log(prob / denominator) )
                        # continue
                        # raise ValueError

                    # else:

                    # beta_range = self.all_beta_ranges[index]
                    # beta_index = np.argmin(np.abs(beta - beta_range))
                    # print(self.cache_denominators[index][beta_index], beta, beta_index, np.max(self.cache_denominators[index]))
                    
                    
                    
                    # if not np.isnan((numerator/denominator)[beta_index]):
                    # print(prob, -np.log(prob))
                    total_ll += -np.log(prob)# *
                             #np.abs(np.diff(self.all_beta_ranges[index]))))
                    # if total_ll == np.inf:
                    #     print(beta, np.sum(llhood), np.sum(numerator), np.sum(denominator),
                    #         self.cache_denominators[index][1:], np.diff(self.all_beta_ranges[index]),
                    #         self.all_weights[index], q_index, self.likelihoods[sample_size][context][q_index, :])
                    #     raise ValueError
                    # print(total_ll)
                    # raise ValueError
                    # total_ll += prob / denominator  # + 1e-8)
                    # print(beta, sd, q_index, prob, denominator, total_ll)
                    # results.append(prob/denominator)

                    # print((prob / denominator + 1e-8), qs, sample_size, closest_sample_size, context, closest_sample_size_inside, closest_context)
        if separate_cs:
            # print(c_positive, c_negative, neutral_prob, total_ll)
            self.param_history.append([c_positive, c_negative, neutral_prob])
            self.likelihood_history.append(total_ll)
        else:
            # print(c, neutral_prob, total_ll)
            self.param_history.append([c, neutral_prob])
            self.likelihood_history.append(total_ll)
        return total_ll#, results

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
