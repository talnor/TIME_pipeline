#!/usr/bin/env python
"""Plot ETI values for calculations based on different fractions of the region of interest in the pol gene.
Used to evaluate the effect of missing data in parts of the region of interest,
i.e. missing data randomly throughout the pol gene, missing data in the beginning or end of the gene, as compared to
the full region"""

import sys
import matplotlib.pyplot as plt
import pandas as pd

file_org = sys.argv[1]
file_random = sys.argv[2]
file_start_missing = sys.argv[3]
file_end_missing = sys.argv[4]
pct = str(sys.argv[5])
anonymize = sys.argv[6]

def get_failed_samples(etifile):
    """Returns a list of samples that lack ETI estimates"""
    df = pd.read_csv(etifile, sep=",", index_col="Sample")
    failed = df[df['ETI'] == "-"]
    samples = failed.index.values.tolist()
    return samples

def get_sample_names(etifile, anonym, outfilename):
    """Returns a list of sample names for use in legends. If legend names are anonymized also output a key file"""
    df = pd.read_csv(etifile, sep=",", index_col="Sample")
    data = df[df['ETI'] != "-"]
    samples = data.index.values.tolist()
    if anonym:
        ids = []
        with open(outfilename, "w") as out:
            for i in range(len(samples)):
                id = f"S{i+1}"
                ids.append(id)
                out.write(f"{samples[i]},{id}\n")
        samples = ids
    return samples

def get_marker(num):
    """Returns a marker type for plotting"""
    marker = "."
    if (num > 6) and (num < 14):
        marker = "x"
    elif (num >= 14) and (num < 21):
        marker = "s"
    elif num >= 21:
        marker = "^"
    return marker + "-"

if pct == "50":
    positions = ["100%", f"random {pct}%", f"2087-3590 ({pct}%)", f"3593-5096 ({pct}%)"]
elif pct == "75":
    positions = ["100%", f"random {pct}%", f"2087-4343 ({pct}%)", f"2840-5096 ({pct}%)"]
else:
    sys.exit("Please enter % of sites used [50, 75].")
files = [file_org, file_random, file_end_missing, file_start_missing]
eti = pd.DataFrame()
failed_samples = get_failed_samples(file_org)
samples = get_sample_names(file_org, anonymize, f"eti_comparison_{pct}pct_region_covered_sampleIDs.txt")
for i in range(len(files)):
    df = pd.read_csv(files[i], sep=",", index_col = "Sample")
    df_filtered = df.drop(index=failed_samples)
    df_filtered = df_filtered.replace("-", -1)
    formatted = df_filtered.astype({"ETI":float}, copy=True, errors="raise")
    y = formatted["ETI"]
    eti[positions[i]] = y

# Create plot
plt.style.use("ggplot")
fig, ax = plt.subplots()
for i in range(len(eti)):
    mark = get_marker(i)
    ax.plot(positions, eti.iloc[i], mark, linewidth=1.0, label=samples[i])
ax.set_ylabel("ETI [years]")
ax.set_xlabel("Percentage of region covered", labelpad=10, fontsize=10)
ax.set_title(f"ETI for subpopulations of {pct}% of sites", fontsize=13, pad=8)
ax.set_xlabel("Part of region covered", labelpad=10, fontsize=10)
ax.legend(bbox_to_anchor=(1.05, 1.0), fontsize="6", loc='upper left')
plt.tight_layout()
plt.savefig(f"eti_comparison_{pct}pct_region_covered.png", dpi=500)