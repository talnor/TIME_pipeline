#!/usr/bin/env python
"""Analyse the shiver BaseFreq files to calculate estimated time of infection using intrasample SNP variation.
Performs the same calculations as the standard ETI script, but does so while specifying a subset of the region of
interest to be used in calculations. What subset of positions should be included is defined by the following parameters:
'mode': ["standard", "random", "missing-beginning", "missing-end"]
'pct_sites_extracted': [0.5, 0.75, 1]

Calculations are based on equation 2 in the following publication:
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
import random

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

# For threshold evaluations
pct_sites_extracted = [0.5, 0.75, 1]
extraction_modes = ["standard", "random", "missing-beginning", "missing-end"]


def get_report_header(param):
    """Return ETI report header"""
    header = ",".join(
        [
            "Project",
            "Sample",
            "Distance",
            "ETI",
            str(param["min_cov"]) + "X_Coverage",
            "1000X_coverage",
            "Average_coverage\n",
        ]
    )
    return header


def format_output_string(
    ticket,
    sample,
    cov_threshold=0,
    cov_1000=0,
    cov_average=0,
    pairwise_distance="-",
    eti="-",
):
    """Return a formatted output string for the results of a sample"""
    sampleline = ",".join(
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
    return sampleline


def create_calculations_report(data, directory, sample, filtered=True):
    """Creates a csv file for the data and calculations performed for a sample"""
    if filtered:
        sample_outfile = os.path.join(
            directory, sample + "_full_filtered_calculations.csv"
        )
    else:
        sample_outfile = os.path.join(
            directory, sample + "_failed_unfiltered_calculations.csv"
        )
    data.to_csv(sample_outfile, sep=",", index=False)


def get_positions_of_interest(param, pct_sites, mode):
    """Return positions in HXB2 coordinates to get data for"""
    positions_of_interest = range(
        param["pos_start"] + param["start_base"] - 1,
        param["pos_fin"] + 1,
        param["step"],
    )
    if mode == "standard":
        positions = positions_of_interest
    # Extract a subset of positions of interest
    elif mode == "random":
        number_of_positions = round(len(positions_of_interest) * pct_sites)
        positions = random.sample(positions_of_interest, number_of_positions)
    elif mode == "missing-beginning":
        positions = positions_of_interest[
            int(len(positions_of_interest) - len(positions_of_interest) * pct_sites) :
        ]
    elif mode == "missing-end":
        positions = positions_of_interest[: int(len(positions_of_interest) * pct_sites)]
    else:
        sys.exit("Error in selection of genome positions")
    return positions


def get_region_data(base_frequency_file, positions):
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
            "N",
        ]
        columns = {
            "HXB2_position": str,
            "Reference_position": str,
            "Reference_Base": str,
        }
        all_data = pd.read_csv(f, header=0, names=headers)
    converted_data = all_data.astype(columns, copy=True, errors="raise")
    pol = converted_data[
        converted_data["HXB2_position"].isin(str(i) for i in positions)
    ]
    return pol


def calculate_nucleotide_frequencies(pol):
    """Normalize nucleotide counts"""
    nucleotide_counts_per_position = pol[["A", "C", "G", "T"]]
    total_counts_per_position = nucleotide_counts_per_position.sum(axis=1)
    pol["total_bases"] = total_counts_per_position.values
    nucleotide_frequency_per_position = nucleotide_counts_per_position.div(
        total_counts_per_position, axis="index"
    )
    nucleotide_frequency_per_position = nucleotide_frequency_per_position.replace(
        [np.inf, -np.inf, "", np.nan], 0
    )
    nucleotide_frequency_per_position.columns = ["freq_A", "freq_C", "freq_G", "freq_T"]
    return pd.concat([pol, nucleotide_frequency_per_position], axis=1)


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
    nucleotide_frequencies_per_position = pol[["freq_A", "freq_C", "freq_G", "freq_T"]]
    diversity = nucleotide_frequencies_per_position * (
        1 - nucleotide_frequencies_per_position
    )
    pol["diversity"] = diversity.sum(axis=1)
    return pol


def calculate_positional_distance(pol, diversity_threshold):
    """Calculate the per base distance"""
    pol = calculate_theta(pol, diversity_threshold)
    pol = calculate_diversity(pol)
    per_base_distance = pol["theta"] * pol["diversity"]
    pol["distance"] = per_base_distance.values
    return pol


def calculate_eti(pairwise_distance, parameters):
    """Calculate estimated time of infection (ETI)"""
    eti = parameters["eti_m"] * pairwise_distance + parameters["eti_c"]
    return eti


def run_analysis(
    samples, eti_summary, outputdir, inputdir, ticket, parameters, pct_sites, mode
):
    """Perform all calculations for a batch of samples and output reports"""
    with open(eti_summary, "w") as out:
        header = get_report_header(parameters)
        out.write(header)
        for sample_info_string in samples:
            sample_info = sample_info_string.strip().lstrip("[").strip("]").strip(",")
            sample = sample_info
            base_frequency_file = glob.glob(
                os.path.join(
                    inputdir, "{}_remap_BaseFreqs_WithHXB2.csv".format(sample_info)
                )
            )
            if len(base_frequency_file) == 0:
                null_data = format_output_string(ticket, sample)
                out.write(null_data)
                continue
            else:
                # Calculate coverage
                positions = get_positions_of_interest(parameters, pct_sites, mode)
                pol = get_region_data(base_frequency_file[0], positions)
                pol = calculate_nucleotide_frequencies(pol)
                cov_average = calculate_average_coverage(pol)
                # Evaluate thresholds
                cov_threshold = check_coverage_threshold(pol, parameters["min_cov"])
                cov_1000 = check_coverage_threshold(pol, 1000)
                pol_unfiltered = pol.copy(deep=True)
                pol = remove_low_coverage_positions(pol, parameters["min_cov"])
                if pol.empty:
                    sample_results = format_output_string(
                        ticket, sample, cov_threshold, cov_1000, cov_average
                    )
                    out.write(sample_results)
                    create_calculations_report(pol_unfiltered, outputdir, sample, False)
                    continue
                else:
                    # Perform ETI calculations
                    pol = calculate_positional_distance(
                        pol, parameters["diversity_threshold"]
                    )
                    pairwise_distance = pol["distance"].mean()
                    pol["avg_pairwise_distance"] = pairwise_distance
                    eti = calculate_eti(pairwise_distance, parameters)
                    pol["ETI"] = eti
                    # Create reports
                    sample_results = format_output_string(
                        ticket,
                        sample,
                        cov_threshold,
                        cov_1000,
                        cov_average,
                        pairwise_distance,
                        eti,
                    )
                    out.write(sample_results)
                    create_calculations_report(pol, outputdir, sample)


def main():
    """Perform eti calculation for different settings for a batch of samples"""
    parameters = {
        "diversity_threshold": diversity_threshold,
        "min_cov": min_cov,
        "pos_start": pos_start,
        "pos_fin": pos_fin,
        "step": step,
        "start_base": start_base,
        "eti_m": eti_m,
        "eti_c": eti_c,
    }

    for mode in extraction_modes:
        for pct_sites in pct_sites_extracted:
            report = "_".join(
                [
                    eti_summary[: eti_summary.rfind(".")],
                    str(mode),
                    "{}PCT.csv".format(int(pct_sites * 100)),
                ]
            )
            outputdir = os.path.join(
                os.path.dirname(eti_summary), "-".join([mode, str(pct_sites * 100)])
            )
            os.mkdir(outputdir)
            run_analysis(
                samples,
                report,
                outputdir,
                inputdir,
                ticket,
                parameters,
                pct_sites,
                mode,
            )


main()
