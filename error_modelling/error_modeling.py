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
Adaptive shrinkage on NaP-TRAP Data.
This code determines a precision-weighted mean and standard error from 4 repeated
experiments of NaP-TRAP. This is then passed to the ASH program to perform
shrinkage (see the "run_ash.r" script).
"""


import pandas as pd
import numpy as np
import scipy.stats as stats


class ErrorModeling:
    def __init__(self):
        self.data = None

    def load(
        self,
        naptrap_results_path,
        predicted_results_path,
        shet_path,
        gene_names_path,
        gnomad_freq_path,
        null_or_nan=False,
        normalization="rpm"
    ):
        """
        Function to load all the prerequisite data

        null_or_nan: if True, add all variants without a predicted effect (i.e. nan) to the Null category.
            If False, exclude these variants
        """
        if normalization != "rpm" and normalization != "spike":
            raise ValueError("normalization must be RPM or spike")
        self.normalization = normalization

        self.null_or_nan = null_or_nan

        # Load csv of input & pulldown for all reporters, including all variants tested
        naptrap_counts = pd.read_csv(naptrap_results_path, index_col=0)

        for replicate in range(1, 5):
            if self.normalization == "rpm":
                # Calculate RPM values
                naptrap_counts.loc[:, "input_rpm_12h_B" + str(replicate)] = (
                    naptrap_counts["input_12h_B" + str(replicate)]
                    / np.sum(naptrap_counts["input_12h_B" + str(replicate)])
                ) * 1000000
                naptrap_counts.loc[:, "pulldown_rpm_12h_B" + str(replicate)] = (
                    naptrap_counts["pulldown_12h_B" + str(replicate)]
                    / np.sum(naptrap_counts["pulldown_12h_B" + str(replicate)])
                ) * 1000000

                # Calculate translation values for both replicates
                naptrap_counts.loc[:, "trans_B" + str(replicate)] = (
                    naptrap_counts["pulldown_rpm_12h_B" + str(replicate)]
                    / naptrap_counts["input_rpm_12h_B" + str(replicate)]
                ) + 0.001
            elif self.normalization == "spike":
                normalizations = naptrap_counts[naptrap_counts.index.str.contains("spk")].sum(axis=0)
                nonspike = ~naptrap_counts.index.str.contains("spk")
                naptrap_counts.loc[:, "input_spike_12h_B" + str(replicate)] = np.full(naptrap_counts.shape[0], np.nan)
                naptrap_counts.loc[:, "pulldown_spike_12h_B" + str(replicate)] = np.full(naptrap_counts.shape[0], np.nan)
                naptrap_counts.loc[:, "input_spike_12h_B" + str(replicate)][nonspike] = (
                    naptrap_counts["input_12h_B" + str(replicate)][nonspike]
                    / normalizations["input_12h_B" + str(replicate)])
                naptrap_counts.loc[:, "pulldown_spike_12h_B" + str(replicate)][nonspike] = (
                    naptrap_counts["pulldown_12h_B" + str(replicate)][nonspike]
                    / normalizations["pulldown_12h_B" + str(replicate)])

                # Calculate translation values for both replicates
                naptrap_counts.loc[:, "trans_B" + str(replicate)] = (
                    naptrap_counts["pulldown_spike_12h_B" + str(replicate)]
                    / naptrap_counts["input_spike_12h_B" + str(replicate)]
                )

        # We'll only use variants with input reads > 0
        naptrap_counts = naptrap_counts[
            (naptrap_counts["input_12h_B1"] > 0)
            & (naptrap_counts["input_12h_B2"] > 0)
            & (naptrap_counts["input_12h_B3"] > 0)
            & (naptrap_counts["input_12h_B4"] > 0)
        ]
        print("After filtering, we have: {} variants".format(naptrap_counts.shape[0]))

        # Load csv of predicted effects for each reporter passing QC
        naptrap_predicted_effects = pd.read_csv(predicted_results_path)

        # Load shet values for each gene
        shet = pd.read_csv(shet_path, delimiter="\t")

        # Load file connecting gene names and hgnc files
        gene_to_hgnc = pd.read_csv(gene_names_path, skiprows=1)

        # Gnomad frequencies for each variant, from Ethan
        gnomad_freq = pd.read_csv(gnomad_freq_path, sep="\t", index_col=0)
        gnomad_freq = gnomad_freq.drop(columns=["Unnamed: 9"])  # drop extraneous column

        # Merge naptrap counts and predicted effects dataframes
        naptrap_predicted_effects["humvar"] = (
            naptrap_predicted_effects.gene_id
            + "_"
            + naptrap_predicted_effects.Chr.astype(str)
            + "-"
            + naptrap_predicted_effects.Position.astype(str)
            + "-"
            + naptrap_predicted_effects.Ref
            + "-"
            + naptrap_predicted_effects.Alt
        )
        naptrap_counts_predicted = pd.merge(
            naptrap_predicted_effects,
            naptrap_counts,
            right_index=True,
            left_on="humvar",
            how="inner",
        )
        out = naptrap_counts[~naptrap_counts.index.isin(naptrap_counts_predicted["humvar"])]
        print(out[~out.index.str.contains('ref')].shape)
        print(naptrap_counts.shape)
        print("After merging with predicted effects we have {} variants".format(naptrap_counts_predicted.shape[0]))
        print("Unique rows {}".format(np.unique(naptrap_counts_predicted["humvar"]).shape[0]))

        # naptrap_counts_predicted = naptrap_counts_predicted.drop(
        #     columns=["ref_seq", "insert_seq"]
        # )

        # Extract gene_id from humvar
        naptrap_counts_predicted["gene_id"] = naptrap_counts_predicted.humvar.str.split(
            "_"
        ).str[0]

        gnomad_freq["gene"] = gnomad_freq["gene"].str.replace("\xa0", "")
        gnomad_freq["gene"] = gnomad_freq["gene"].str.replace(".00", "")
        gnomad_freq["gene"] = gnomad_freq["gene"].str.replace(".00", "")
        gnomad_freq["humvar"] = gnomad_freq.gene + "_" + gnomad_freq.variant_name

        merged = pd.merge(
            gnomad_freq,
            naptrap_counts_predicted,
            on=["humvar"], how="inner"
        )
        print("After merging with frequencies we have {} variants".format(merged.shape[0]))
        print("Unique rows {}".format(np.unique(merged["humvar"]).shape[0]))

        merged = pd.merge(
            gene_to_hgnc, merged, left_on=["Input"], right_on="gene", how="inner"
        )
        print("After merging with hgnc we have {} variants".format(merged.shape[0]))
        print("Unique rows {}".format(np.unique(merged["humvar"]).shape[0]))
        assert np.sum(np.isnan(merged["allele_count"]) == 0)

        self.data = pd.merge(merged, shet, left_on="HGNC ID", right_on="hgnc", how="inner")
        print("After merging with shet we have {} variants".format(self.data.shape[0]))
        print("Unique rows {}".format(np.unique(self.data["humvar"]).shape[0]))

    def precision_weighted_mean_se(self, output_path):
        # Determine the median translation values for each gene
        median_norm_1 = np.full(self.data.shape[0], np.nan, dtype=float)
        median_norm_2 = np.full(self.data.shape[0], np.nan, dtype=float)
        median_norm_3 = np.full(self.data.shape[0], np.nan, dtype=float)
        median_norm_4 = np.full(self.data.shape[0], np.nan, dtype=float)

        # Remove genes with a median translation value of 0
        removed_genes = []
        # "Input" is the gene name
        for gene in np.unique(self.data.Input):
            gene_vars = self.data[self.data.Input == gene]
            median_norm_1[self.data.Input == gene] = gene_vars["trans_B1"] / np.median(
                gene_vars["trans_B1"]
            )
            median_norm_2[self.data.Input == gene] = gene_vars["trans_B2"] / np.median(
                gene_vars["trans_B2"]
            )
            median_norm_3[self.data.Input == gene] = gene_vars["trans_B3"] / np.median(
                gene_vars["trans_B3"]
            )
            median_norm_4[self.data.Input == gene] = gene_vars["trans_B4"] / np.median(
                gene_vars["trans_B4"]
            )

            if np.any(
                np.array(
                    [
                        np.median(gene_vars["trans_B1"]),
                        np.median(gene_vars["trans_B2"]),
                        np.median(gene_vars["trans_B3"]),
                        np.median(gene_vars["trans_B4"]),
                    ]
                )
                == 0
            ):
                removed_genes.append(gene)

        for replicate, median_norm in enumerate(
            [median_norm_1, median_norm_2, median_norm_3, median_norm_4]
        ):
            # Haploid linear scale
            self.data.loc[:, "norm_B" + str(replicate + 1)] = median_norm - 1

            # # Convert to diploid scale
            # self.data.loc[:, "norm_B" + str(replicate + 1)] = np.log2(
            #     1 + self.data["norm_B" + str(replicate + 1)] / 2
            # )

        # Drop genes with median translation values of 0
        print("Before dropping median values of 0", self.data.shape)
        self.data = self.data[~np.isin(self.data.gene, removed_genes)]
        print("After dropping median values of 0", self.data.shape)

        columns_and_new_names = {
            "trans_B1": "expected_pulldown_B1",
            "trans_B2": "expected_pulldown_B2",
            "trans_B3": "expected_pulldown_B3",
            "trans_B4": "expected_pulldown_B4",
        }

        # Iterate over each original column and its corresponding new column name
        for replicate, (orig_col, new_col) in enumerate(columns_and_new_names.items()):
            # Group by 'gene' and compute the median for each group
            self.data[new_col] = (
                self.data.groupby("gene")[orig_col].transform("median")
                * self.data["input_" + self.normalization + "_12h_B" + str(replicate + 1)]
            )

        # for gene in np.unique(self.data.gene_id):
        #     gene_vars = self.data[self.data.gene_id == gene]
        #     self.data[self.data.gene_id == gene, "gene_trans_B" + str(replicate) + "_median"] = np.median(gene_vars["trans_B" + str(replicate_inner)])
        #     for replicate_inner in range(1, 5):
        #         expected[self.data.gene_id == gene] = (
        #             np.median(gene_vars["trans_B" + str(replicate_inner)])
        #             * gene_vars["input_12h_B" + str(replicate_inner)]
        #         )

        # for replicate_outer in range(1, 5):
        #     # Calculate the median translation value for each gene and compute expected pulldown
        #     expected = np.full(self.data.shape[0], np.nan, dtype=float)
        # self.data.loc[:, "expected_pulldown_B" + str(replicate_outer)] = expected

       
        # Bin the genes based on expected pulldown values (using replicate 1 as an example, adjust as necessary)
        probs = np.linspace(0.01, 0.99, 100)

        for replicate in range(1, 5):
            bin_values = stats.mstats.mquantiles(
                self.data["expected_pulldown_B" + str(replicate)], prob=probs
            )

            bin_values = (
                [self.data["expected_pulldown_B" + str(replicate)].min()]
                + list(bin_values)
                + [self.data["expected_pulldown_B" + str(replicate)].max() + 1e-10]
            )

            self.data.loc[:, "bins_B" + str(replicate)] = np.digitize(
                self.data["expected_pulldown_B" + str(replicate)], bin_values
            )
            SE_arr = np.full(self.data.shape[0], np.nan, dtype=float)
            for bin in range(1, len(bin_values)):
                subset = self.data[self.data["bins_B" + str(replicate)] == bin]
                # Find all 4 choose 2 covariances and take the mean
                cov_1_2 = np.cov(subset["norm_B1"], subset["norm_B2"], ddof=1)[0][1]
                cov_1_3 = np.cov(subset["norm_B1"], subset["norm_B3"], ddof=1)[0][1]
                cov_1_4 = np.cov(subset["norm_B1"], subset["norm_B4"], ddof=1)[0][1]
                cov_2_3 = np.cov(subset["norm_B2"], subset["norm_B3"], ddof=1)[0][1]
                cov_2_4 = np.cov(subset["norm_B2"], subset["norm_B4"], ddof=1)[0][1]
                cov_3_4 = np.cov(subset["norm_B3"], subset["norm_B4"], ddof=1)[0][1]

                overall_cov = np.mean(
                    [cov_1_2, cov_1_3, cov_1_4, cov_2_3, cov_2_4, cov_3_4]
                )
                # SE = np.sqrt((np.var(subset["norm_B" + str(replicate)]) ) / subset.shape[0])
                # Var = np.sqrt(
                #     np.sum(
                #         (
                #             subset["norm_B" + str(replicate)]
                #             - np.mean(subset["norm_B" + str(replicate)])
                #         )
                #         ** 2
                #     )
                #     / subset.shape[0]
                # )
                SE = np.sum(
                        (
                            subset["norm_B" + str(replicate)]
                            - np.mean(subset["norm_B" + str(replicate)])
                        )
                        ** 2
                    ) / subset.shape[0]

                SE = SE - overall_cov

                # SE = np.sqrt(SE/subset.shape[0])        
            
                print(f"Bin {bin} Variance for trans_b {replicate}: {SE}")
                SE_arr[self.data["bins_B" + str(replicate)] == bin] = SE

                # SE_arr[self.data["bins_B" + str(replicate)] == bin] = (
                #     SE - self.overall_cov
                # )

            self.data.loc[:, "SE_B" + str(replicate)] = SE_arr

            print(
                "Variance of B"
                + str(replicate)
                + ": "
                + str(self.data["norm_B" + str(replicate)].var())
            )

        # precision weighted mean
        num = (
            (self.data["norm_B1"] / self.data["SE_B1"] ** 2)
            + (self.data["norm_B2"] / self.data["SE_B2"] ** 2)
            + (self.data["norm_B3"] / self.data["SE_B3"] ** 2)
            + (self.data["norm_B4"] / self.data["SE_B4"] ** 2)
        )
        denom = (
            self.data["SE_B1"] ** (-2)
            + self.data["SE_B2"] ** (-2)
            + self.data["SE_B3"] ** (-2)
            + self.data["SE_B4"] ** (-2)
        )
        precision_weighted_mean = num / denom
        precision_weighted_se = np.sqrt(1 / denom)
        self.data["prec_weighted_mean"] = precision_weighted_mean
        self.data["prec_weighted_se"] = precision_weighted_se

        # All variants
        pd.DataFrame(
            {"mean": precision_weighted_mean, "ses": precision_weighted_se}
        ).to_csv(output_path + "_all_predicted.csv")

        # Negative predicted variants
        pd.DataFrame(
            {
                "mean": precision_weighted_mean[
                    self.data["predicted_category"] == "Negative"
                ],
                "ses": precision_weighted_se[
                    self.data["predicted_category"] == "Negative"
                ],
            }
        ).to_csv(output_path + "_negative_predicted.csv")

        # Null predicted variants
        if self.null_or_nan:
            pd.DataFrame(
                {
                    "mean": precision_weighted_mean[
                        np.logical_or(
                            self.data["predicted_category"] == "Null",
                            self.data["predicted_category"].astype(str) == "nan",
                        )
                    ],
                    "ses": precision_weighted_se[
                        np.logical_or(
                            self.data["predicted_category"] == "Null",
                            self.data["predicted_category"].astype(str) == "nan",
                        )
                    ],
                }
            ).to_csv(output_path + "_null_predicted.csv")
        else:
            pd.DataFrame(
                {
                    "mean": precision_weighted_mean[
                        self.data["predicted_category"] == "Null"
                    ],
                    "ses": precision_weighted_se[
                        self.data["predicted_category"] == "Null"
                    ],
                }
            ).to_csv(output_path + "_null_predicted.csv")

        # Positive predicted variants
        pd.DataFrame(
            {
                "mean": precision_weighted_mean[
                    self.data["predicted_category"] == "Positive"
                ],
                "ses": precision_weighted_se[
                    self.data["predicted_category"] == "Positive"
                ],
            }
        ).to_csv(output_path + "_positive_predicted.csv")

        # All variants
        pd.DataFrame(
            {
                "mean": precision_weighted_mean,
                "ses": precision_weighted_se,
            }
        ).to_csv(output_path + "_all_predicted.csv")
