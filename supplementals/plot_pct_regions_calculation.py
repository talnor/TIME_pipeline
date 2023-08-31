#!/usr/bin/env python
"""Plot the ETI values for calculations based on 100%, 75%, and 50% of the positions in the region of interest.
Specify the mode of selection of sites performed during calculations as the
'mode': [random, missing-beginning, missing-end]"""

import sys
import matplotlib.pyplot as plt
import pandas as pd

file_org = sys.argv[1]
file_50 = sys.argv[2]
file_75 = sys.argv[3]
mode = sys.argv[4]
anonymize = sys.argv[5]


def get_failed_samples(etifile):
    """Returns a list of samples that lack ETI estimates"""
    df = pd.read_csv(etifile, sep=",", index_col="Sample")
    failed = df[df["ETI"] == "-"]
    samples = failed.index.values.tolist()
    return samples


def get_sample_names(etifile, anonym, outfilename):
    """Returns a list of sample names for use in legends. If legend names are anonymized also output a key file"""
    df = pd.read_csv(etifile, sep=",", index_col="Sample")
    data = df[df["ETI"] != "-"]
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


outfile = f"eti_comparison_{mode}_region_covered.png"
if mode == "random":
    title = "ETI for random positions"
    positions = ["100%", f"{mode} 75%", f"{mode} 50%"]
elif mode == "missing-beginning":
    title = "ETI for last stretch of positions"
    mode = "last"
    positions = ["100%", "2840-5096 (75%)", "3593-5096 (50%)"]
elif mode == "missing-end":
    title = "ETI for first stretch of positions"
    mode = "first"
    positions = ["100%", "2087-4343 (75%)", "2087-3590 (50%)"]
else:
    sys.exit("Please enter mode type [random, missing-beginning, missing-end].")

files = [file_org, file_75, file_50]
eti = pd.DataFrame()
failed_samples = get_failed_samples(file_org)
samples = get_sample_names(
    file_org, anonymize, f"eti_comparison_{mode}_region_covered_sampleIDs.txt"
)
for i in range(len(files)):
    df = pd.read_csv(files[i], sep=",", index_col="Sample")
    df = df.drop(index=failed_samples)
    df = df.replace("-", -1)
    formatted = df.astype({"ETI": float}, copy=True, errors="raise")
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
ax.set_title(title, fontsize=13, pad=8)
ax.set_xlabel("Part of region covered", labelpad=10, fontsize=10)
ax.legend(bbox_to_anchor=(1.05, 1.0), fontsize="6", loc="upper left")
plt.tight_layout()
plt.savefig(outfile, dpi=500)
