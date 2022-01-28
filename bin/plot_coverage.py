#!/usr/bin/env python
"""Plot the genome coverage"""

import sys
import matplotlib.pyplot as plt
import pandas as pd

# Get data
inputf = sys.argv[1]
sample = sys.argv[2]
df = pd.read_csv(inputf, sep='\t')

# Plot data
fig, ax = plt.subplots()
ax.plot(df.iloc[:, 1], df.iloc[:, 2])

# Add descriptors
ax.set_xlabel('Reference position', labelpad=8)
ax.set_ylabel('Coverage [bp]', labelpad=8)
ax.set_title(sample, fontsize=14, pad=8)

# Output plot
outlog = inputf[:-3] + 'png'
plt.savefig(outlog)

# Also make log-plot
plt.yscale('log')
outlog = inputf[0:-3] + 'log.png'
plt.savefig(outlog)
