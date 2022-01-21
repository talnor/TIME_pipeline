import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn

# Get data
inputf = sys.argv[1]
sample = sys.argv[2]
df = pd.read_csv(inputf, sep='\t')

# Set plot style
plt.style.use('seaborn')
#seaborn.set(rc={"xtick.bottom" : True, "ytick.left" : True})

# Plot data
fig, ax = plt.subplots()
ax.plot(df.iloc[:,1], df.iloc[:,2])

# Add descriptors
ax.set_xlabel('Reference position', labelpad=8)
ax.set_ylabel('Coverage [bp]', labelpad=8)
ax.set_title(sample, fontsize=14, pad=8)

# Output plot
outlog = inputf[0:-3] + 'png'
plt.savefig(outlog)
