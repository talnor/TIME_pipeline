"""Script to calculate mean identity for a batch of samples. Runs two tools by Shiver - LinkIdentityToCoverage.py and
LinkIdentityToCoverage_CombineBams.py in addition to separate plotting of samples.
Both scripts are available at the Shiver github https://github.com/ChrisHIV/shiver/blob/master/tools/
Wymant, C. et al. (2016). Easy and accurate reconstruction of whole hiv genomes from short-read sequence data. bioRxiv.5
Uses a bam file, reference file and a start and stop position in reference coordinates to calculate identity.
"""

import os
import sys
import glob
import pandas as pd
import subprocess
import matplotlib.pyplot as plt

inputdir = sys.argv[1]
anonymize = sys.argv[2]
samples = sys.argv[3:]

start_HXB2 = 2085
stop_HXB2 = 5096

def get_results_filename(inputdir, sample, resultsfile):
    """Return full filename for type of results specified"""
    name_dict = {"bam": "_remap.bam", "ref":"_remap_ref.fasta", "basefreq": "_remap_BaseFreqs_WithHXB2.csv"}
    infile = os.path.join(inputdir, sample, "shiver", sample + name_dict[resultsfile])
    full_path = glob.glob(infile)
    if len(full_path) == 0:
        infile = os.path.join(inputdir, sample, sample + name_dict[resultsfile])
        full_path = glob.glob(infile)
    if len(full_path) == 0:
        infile = os.path.join(inputdir, sample + name_dict[resultsfile])
        full_path = glob.glob(infile)
    return full_path

def get_reference_positions(base_frequency_file, start_HXB2, stop_HXB2):
    """Returns start and stop positions in the coordinates of the custom reference used for mapping"""
    with open(base_frequency_file, "r") as f:
        df = pd.read_csv(f, header=0, index_col=0, usecols=[0, 1])
    start_ref = df.loc[str(start_HXB2)][0]
    stop_ref = df.loc[str(stop_HXB2)][0]
    return start_ref, stop_ref

def plot_identity_per_sample(csvfiles, sampleID, anonymize):
    """Plot the mean read identity for each sample separately
    Input files have the file structure "Coverage,Number of positions with that coverage,Mean identity"
    """
    # Create plot
    plt.style.use("ggplot")
    fig, ax = plt.subplots()
    i = 0
    marker = "."
    for csv in csvfiles:
        # Set marker
        i += 1
        if (i > 7) and (i < 15):
            marker = "x"
        elif (i >= 15) and (i < 22):
            marker = "*"
        elif i >= 22:
            marker = "+"
        # Plot sample data
        sample = csv[:csv.rfind("_")]
        if anonymize:
            sample = sampleID[sample]
        else:
            sample = csv.split("_")[2]
        df = pd.read_csv(csv, sep=",", usecols = ["Coverage", "Mean identity"])
        formatted = df.astype({"Coverage":float, "Mean identity":float}, copy=True, errors="raise")
        ax.plot(formatted["Coverage"], formatted["Mean identity"], marker, linewidth=1.0, label=sample)
    # Format plot
    ax.set_ylabel("Mean identity")
    ax.set_xlabel("Coverage", labelpad=10, fontsize=10)
    ax.set_title("The read identity averaged over all genome positions for each sample", fontsize="10", pad=8)
    ax.legend(bbox_to_anchor=(1.05, 1.0), fontsize="6", loc='upper left')
    plt.tight_layout()
    plt.savefig("per_sample_mean_identity.png", dpi=500)
    plt.xscale("log")
    plt.savefig("per_sample_mean_identity_log.png", dpi=500)
    plt.close()

def plot_identity_samples_averaged(csvfiles, sampleID, anonymize, single_sample=False):
    """Use the Shiver tool LinkIdentityToCoverage_CombineBams.py to plot the mean read identity averaged over all genome positions in all samples supplied"""
    if single_sample:
        if anonymize:
            single_sample = sampleID[single_sample]
        outfile = f"identity_coverage_linkage_{single_sample}"
        plottitle = f"The read identity averaged over all genome positions for sample {single_sample}"
    else:
        outfile = "identity_coverage_linkage_all_samples"
        plottitle = "The read identity averaged over all genome positions in all samples"
    proc = subprocess.run(["LinkIdentityToCoverage_CombineBams.py", "--title", plottitle, '--x-min-max', "1,10000", '--y-min-max', "0.70,1.0", outfile] + csvfiles)

def run_analysis(samples, inputdir, start_HXB2, stop_HXB2, anonymize):
    """Perform all calculations for a batch of samples"""
    csvfiles = []
    sampleID = {}
    if anonymize:
        for i in range(len(samples)):
            sampleID[samples[i]] = f"S{i+1}"
        with open("sampleIDs.txt", "w") as out:
            for sample, ID in sampleID.items():
                out.write(f"{sample},{ID}\n")
    for sample_info_string in samples:
        sample_info = sample_info_string.strip().lstrip("[").strip("]").strip(",")
        sample = sample_info
        base_frequency_file = get_results_filename(inputdir, sample, "basefreq")
        if len(base_frequency_file) != 0:
            start_ref, stop_ref = get_reference_positions(base_frequency_file[0], start_HXB2, stop_HXB2)
            bam_file = get_results_filename(inputdir, sample, "bam")[0]
            ref_file = get_results_filename(inputdir, sample, "ref")[0]
            outfile = f"{sample}_identities.csv"
            csvfiles.append(outfile)
            proc = subprocess.run(["LinkIdentityToCoverage.py", bam_file, ref_file, "--start", str(start_ref), "--end", str(stop_ref)], stdout = subprocess.PIPE, encoding='utf-8')
            with open(outfile, "w") as out:
                out.write(str(proc.stdout))
            plot_identity_samples_averaged([outfile], sampleID, anonymize, sample)
        else:
            continue
    plot_identity_samples_averaged(csvfiles, sampleID, anonymize)
    plot_identity_per_sample(csvfiles, sampleID, anonymize)



run_analysis(samples, inputdir, start_HXB2, stop_HXB2, anonymize)