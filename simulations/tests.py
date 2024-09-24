"""
Test cases for the selection inference algorithm
"""
from tqdm import tqdm
import collections
import json
import logging
import unittest

import core

import numpy as np
import pytest

def sample_q_jeff(s, contexts, sample_sizes, ll_dict, neutral_prop=0, unobserved=False, progress=False):
    # Sample allele frequencies using given likelihoods
    # Works for one mutation rate at a time
    if unobserved == False:
        assert ll_dict.normalized == True

    closest_sample_sizes = []
    qs = []
    probs_list = []
    sampled_probs = 0

    neutral_sites = np.random.choice(np.arange(0, s.shape[0]), int(neutral_prop * s.shape[0]), replace=False)
    print(neutral_sites)

    for index, (s_cur, context, sample_size) in tqdm(enumerate(zip(s, contexts, sample_sizes)), total=s.shape[0], disable=not progress):
        closest_s = np.argmin(np.abs(s_cur - ll_dict.S_GRID))
        sample_size_keys = np.array(list(ll_dict.likelihoods.keys()))
        if sample_size >= np.max(sample_size_keys):
            closest_sample_size = sample_size_keys[-1]
        elif sample_size in sample_size_keys:
            closest_sample_size = sample_size
        else:
            closest_sample_size = sample_size_keys[
                np.searchsorted(
                    sample_size_keys,
                    [
                        sample_size,
                    ],
                    side="right",
                )[0]
            ]
        closest_sample_sizes.append(closest_sample_size)

        if unobserved:
            allele_frequencies = np.arange(0, closest_sample_size + 1)
        else:
            allele_frequencies = np.arange(1, closest_sample_size)
        if index in neutral_sites:
            sampled_q = np.random.choice(allele_frequencies,
                                   p=ll_dict.likelihoods[closest_sample_size][context][:, 0])
        else:
            sampled_q = np.random.choice(allele_frequencies,
                                   p=ll_dict.likelihoods[closest_sample_size][context][:, closest_s])
        # Because this is a normalized dictionary, we need to subtract 1 to get correct index
        #sampled_prob = ll_dict.likelihoods[closest_sample_size][context][sampled_q - 1, closest_s]
        #sampled_probs += sampled_prob
        #probs_list.append(sampled_prob)
        assert sampled_q <= closest_sample_size
        qs.append(sampled_q)
    return np.array(qs), sampled_probs, closest_sample_sizes, probs_list



class TestBasicFunctions:
    """
    Test accuracy of inference
    """

    def create_ll(self):
        S_GRID = np.array(
           [0.] + np.exp(np.linspace(np.log(1e-8), 0, num=100)).tolist()
        )
        ll_test = core.Likelihoods({150000: "lls_subsampled_150k.pkl"}, S_GRID[:-1], "mutation_rate_methylation_bins.txt")
        ll_test.load(normalize=True)
        return ll_test


    def test_basic(self):
        c = 1
        sample_size = 15000
        beta = 0.1
        sd = 0
        s_het = 0.1
        context = "A(C>G)C_0A"
        ll = self.create_ll()

        #assert ll.compute_ll(c, sample_size, beta, sd, s_het, context) == 



class TestPieceWise:
    """
    Test piecewise inference of c
    """
    def create_ll(self):
        S_GRID = np.array(
           [0.] + np.exp(np.linspace(np.log(1e-8), 0, num=100)).tolist()
        )
        ll_test = core.Likelihoods({150000: "lls_subsampled_150k.pkl"}, S_GRID[:-1], "mutation_rate_methylation_bins.txt")
        ll_test.load(normalize=True)
        return ll_test

    def test_no_sd_no_neutral(self):
        ll_test = self.create_ll()
        c = 1
        neutral_prob = 0
        sample_size = [150000, 150000]
        beta = [-0.1, 0.1]
        sd = [0, 0]
        s_het = [ll_test.S_GRID[0], ll_test.S_GRID[1]]
        q = [1, 2]
        context = ["A(C>G)C_0A", "A(C>G)C_0A"]
        ll_test.load_data(sample_size, beta, sd, s_het, q, context)
        ll_val = ll_test.compute_ll_piecewise_bin_no_se_no_neutral(1, -1, 1)
        # hand compute value to test accuracy
        hand_calc_val = -(np.log(ll_test.likelihoods[150000][context[0]][0,0]) + np.log(ll_test.likelihoods[150000][context[0]][1,1]))
        assert ll_val == hand_calc_val
        ll_val = ll_test.compute_ll_piecewise_bin_no_se_no_neutral(1, -1, 0)
        hand_calc_val = -(np.log(ll_test.likelihoods[150000][context[0]][0,0]))
        assert ll_val == hand_calc_val
        ll_val = ll_test.compute_ll_piecewise_bin_no_se_no_neutral(2, 0, 1)
        np.argmin(np.abs(ll_test.S_GRID - (ll_test.S_GRID[1] * 2)))
        hand_calc_val = -(np.log(ll_test.likelihoods[150000][context[0]][1,5]))
        assert ll_val == hand_calc_val


# Todo: test total probability of larger models   
