"""
Test cases for the selection inference algorithm
"""
import collections
import json
import logging
import unittest

import core

import numpy as np
import pytest



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
        ll_val = ll_test.compute_ll_piecewise_bin_no_se(1, -1, 1)
        # hand compute value to test accuracy
        hand_calc_val = -(np.log(ll_test.likelihoods[150000][context[0]][0,0]) + np.log(ll_test.likelihoods[150000][context[0]][1,1]))
        assert ll_val == hand_calc_val
        ll_val = ll_test.compute_ll_piecewise_bin_no_se(1, -1, 0)
        hand_calc_val = -(np.log(ll_test.likelihoods[150000][context[0]][0,0]))
        assert ll_val == hand_calc_val
        ll_val = ll_test.compute_ll_piecewise_bin_no_se(2, 0, 1)
        np.argmin(np.abs(ll_test.S_GRID - (ll_test.S_GRID[1] * 2)))
        hand_calc_val = -(np.log(ll_test.likelihoods[150000][context[0]][1,5]))
        assert ll_val == hand_calc_val

# Todo: test total probability of larger models   
