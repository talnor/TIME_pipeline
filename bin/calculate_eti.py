#!/usr/bin/env python
"""Analyse the shiver BaseFreq files to calculate estimated time of infection using intrapatient SNP variation.

Calculcations are based on equation 2 in the following publication:
Puller V, Neher R, Albert J (2017), Estimating time of HIV-1 infection from next-generation sequence diversity.
PLoS Comput Biol 13(10): e1005775. https://doi.org/10.1371/journal.pcbi.1005775
As well as using ETI calculations from:
https://hiv.biozentrum.unibas.ch/ETI/
"""

import os
import sys
import glob
import pandas as pd
import numpy as np

inputdir = sys.argv[1]
coverage_threshold = sys.argv[2]
ticket = sys.argv[3]
eti_summary = sys.argv[4]
samples = sys.argv[5:]

# diversity threshold (fraction) for distance calculations
diversity_threshold = 0.01

# Per base coverage threshold for inclusion
min_cov = int(coverage_threshold)

# region
pos_start = 2085
pos_fin = 5096

# ETI calculation formula (eti = mx + c)
eti_m = 259
eti_c = 0.05

# Use every n'th base
step = 3

# Start base in codon
start_base = 3


def get_region_data(base_frequency_file, pos_start, start_base, pos_fin, step):
    """Return base count data for region of interest"""
    with open(base_frequency_file, "r") as f:
        headers = [
            "HXB2_position",
            "Reference_position",
            "Reference_Base",
            "A",
            "C",
            "G",
            "T",
            "gap",
            "N"
        ]
        columns = {
            "HXB2_position": str,
            "Reference_position": str,
            "Reference_Base": str,
            "A": float,
            "C": float,
            "G": float,
            "T": float,
            "gap": float,
            "N": float,
        }
        all_data = pd.read_csv(f, header=0, names=headers)
        all_data.astype(columns)
    positions_of_interest = range(pos_start + start_base - 1, pos_fin + 1, step)
    pol = all_data[all_data["HXB2_position"].isin(str(i) for i in positions_of_interest)]
    return pol


def output_template_data(data, ticket, sample, cov_average, cov_threshold=0, cov_1000=0):
    """Output template data for samples with insuffiecient or 0 coverage in entire region of interest"""
    template_data = [
        ticket,
        sample,
        "-",
        "-",
        str(round(cov_threshold, 2)) + "%",
        str(round(cov_1000, 2)) + "%",
        str(int(round(cov_average))) + "\n",
    ]
    out.write(",".join(template_data))
    sample_outfile = sample + "_pairwise_distance_low_coverage.csv"
    data.to_csv(sample_outfile, sep=",", index=False)


def calculate_nucleotide_frequencies(pol):
    """Normalize nucleotide counts"""
    nucleotide_counts_per_base = pol[["A", "C", "G", "T"]]
    total_counts_per_base = nucleotide_counts_per_base.sum(axis=1)
    pol["total_bases"] = total_counts_per_base.values
    nucleotide_frequency_per_base = nucleotide_counts_per_base.div(
        total_counts_per_base, axis="index"
    )
    nucleotide_frequency_per_base = nucleotide_frequency_per_base.replace(
        [np.inf, -np.inf, "", np.nan], 0
    )
    nucleotide_frequency_per_base.columns = ["freq_A", "freq_C", "freq_G", "freq_T"]
    return pd.concat([pol, nucleotide_frequency_per_base], axis=1)


def calculate_average_coverage(pol):
    """Calculate average coverage for a sample"""
    return pol["total_bases"].mean()


def check_coverage_threshold(pol, min_cov):
    """Check if coverage for each position is above threshold.
    Returns percentage of bases that pass the coverage threshold"""
    column = str(min_cov) + "X_Coverage"
    pol.loc[pol["total_bases"] >= min_cov, column] = "Yes"
    pol.loc[pol["total_bases"] < min_cov, column] = "No"
    if "Yes" in pol[column].values:
        pct_coverage = pol[column].value_counts()["Yes"] / len(pol) * 100
    else:
        pct_coverage = 0
    return pct_coverage


def remove_low_coverage_positions(pol, min_cov):
    """Remove all positions with insufficient read depth"""
    pol = pol[pol[str(min_cov) + "X_Coverage"] == "Yes"]
    return pol


def calculate_theta(pol, diversity_threshold):
    """Calculate theta value for each position"""
    major_base_frequency = pol[["freq_A", "freq_C", "freq_G", "freq_T"]].max(axis=1)
    pol["freq_major_base"] = major_base_frequency.values
    theta_array = np.where((1 - pol["freq_major_base"]) > diversity_threshold, 1, 0)
    pol["theta"] = theta_array.tolist()
    pol.loc[pol["freq_major_base"] == 0, "theta"] = 0
    return pol


def calculate_diversity(pol):
    """Calculate nucleotide diversity value for the sample"""
    nucleotide_frequencies_per_base = pol[["freq_A", "freq_C", "freq_G", "freq_T"]]
    diversity = nucleotide_frequencies_per_base * (1 - nucleotide_frequencies_per_base)
    pol["diversity"] = diversity.sum(axis=1)
    return pol


def calculate_positional_distance(pol, diversity_threshold):
    """Calculate the per base distance"""
    pol = calculate_theta(pol, diversity_threshold)
    pol = calculate_diversity(pol)
    per_base_distance = pol["theta"] * pol["diversity"]
    pol["distance"] = per_base_distance.values
    return pol


def calculate_ETI(pairwise_distance, eti_m, eti_c):
    """Calculate estimated time of infection (ETI)"""
    eti = eti_m * pairwise_distance + eti_c
    return eti


with open(eti_summary, "w") as out:
    out.write(
        ",".join(
            [
                "Project",
                "Sample",
                "Distance",
                "ETI",
                str(min_cov) + "X_Coverage",
                "1000X_coverage",
                "Average_coverage\n",
            ]
        )
    )
    for sample_info in samples:
        sample_info = sample_info.strip().lstrip("[").strip("]").strip(",")
        sample = sample_info.split("_")[1]
        base_frequency_file = glob.glob(
            os.path.join(
                inputdir, "{}_remap_BaseFreqs_WithHXB2.csv".format(sample_info)
            )
        )
        if len(base_frequency_file) == 0:
            template_data = [ticket, sample, "-", "-", "0%", "0%", "0\n"]
            out.write(",".join(template_data))
            continue
        else:
            pol = get_region_data(
                base_frequency_file[0], pos_start, start_base, pos_fin, step
            )
            pol = calculate_nucleotide_frequencies(pol)
            cov_average = calculate_average_coverage(pol)
            cov_threshold = check_coverage_threshold(pol, min_cov)
            cov_1000 = check_coverage_threshold(pol, 1000)
            pol_unfiltered = pol.copy(deep=True)
            pol = remove_low_coverage_positions(pol, min_cov)
            if pol.empty:
                output_template_data(pol_unfiltered, ticket, sample, cov_average, cov_threshold, cov_1000)
                continue
            else:
                pol = calculate_positional_distance(pol, diversity_threshold)
                pairwise_distance = pol["distance"].mean()
                pol["avg_pairwise_distance"] = pairwise_distance
                eti = calculate_ETI(pairwise_distance, eti_m, eti_c)
                pol["ETI"] = eti
                out.write(
                    ",".join(
                        [
                            ticket,
                            sample,
                            str(pairwise_distance),
                            str(eti),
                            str(round(cov_threshold, 2)) + "%",
                            str(round(cov_1000, 2)) + "%",
                            str(int(round(cov_average))) + "\n",
                        ]
                    )
                )
                sample_outfile = sample + "_pairwise_distance.csv"
                pol.to_csv(sample_outfile, sep=",", index=False)
