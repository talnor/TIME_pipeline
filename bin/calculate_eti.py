"""Analyse the shiver BaseFreq files for 3 base variation to get estimates of intrapatient SNP variation"""

import os
import sys
import glob
import pandas as pd
import numpy as np

inputdir = sys.argv[1]	#folder for the project
coverage_threshold = sys.argv[2]

# diversity threshold (fraction) for distance calculations
tc = 0.01

# Per base coverage threshold for inclusion
min_cov = coverage_threshold

# region
pos_start = 2085
pos_fin = 5096

# ETI calculation formula (eti = mx + c)
eti_m = 259
eti_c = 0.05

# Use every n'th base
n = 3

# Start base in codon
base = 3

#Get all the subdirectories i.e. sample directories
subdirs = []
for cont in os.listdir(inputdir):
    if os.path.isdir(os.path.join(inputdir, cont)):
        subdirs.append(cont) 

# Analyse each sample in the directory
print('Project', 'Sample', 'Date', 'Distance', 'ETI', str(min_cov) + 'X_Coverage', '1000X_coverage', 'Average_coverage', sep=',')
for sampledir in subdirs:

    # Get the project name
    proj_name = inputdir.strip().strip('/').split('/')[-1].split('_')[0]

    # Get the sample base frequency file
    files = sorted(glob.glob(os.path.join(inputdir, sampledir, 'shiver*', "*_remap_BaseFreqs_WithHXB2.csv")))

    # No results found    
    if not files:
        print(proj_name, sampledir, *(['-']*6), sep=',')
        continue

    # Calculate ETI
    else:     
        # Get the base of the file name
        SID = os.path.basename(files[0])[:-4]
        outfile = os.path.join(os.path.dirname(files[0]), SID + '_pairwise_distance.csv')
        
        with open(files[0], 'r') as f:
 
            # read the csv    
            df = pd.read_csv(f)
            df.columns = ['HXB2_position', 'Reference_position', 'Reference_Base', 'A', 'C', 'G', 'T', 'gap', 'N']

            # Extract the third nucleotide for each codon in the pol gene region
            pol = df[df['HXB2_position'].isin(str(i) for i in range(pos_start + base - 1, pos_fin + 1, n))]

            # Normalize nucleotide counts
            pol_nucl = pol[['A', 'C', 'G', 'T']]
            tot = pol_nucl.sum(axis=1)
            pol_nucl = pol_nucl.div(tot, axis='index')
            pol_nucl = pol_nucl.replace([np.inf, -np.inf, '', np.nan], 0)
            pol_nucl.columns = ['freq_A', 'freq_C', 'freq_G', 'freq_T']
            pol = pd.concat([pol, pol_nucl], axis=1)
            
            # Check coverage
            avg_cov = tot.mean()
            pol.loc[tot >= min_cov, str(min_cov) + 'X_Coverage'] = 'Yes'
            pol.loc[tot < min_cov, str(min_cov) + 'X_Coverage'] = 'No'
            pol.loc[tot >= 1000, '1000X_Coverage'] = 'Yes'
            pol.loc[tot < 1000, '1000X_Coverage'] = 'No'
            if 'Yes' in pol[str(min_cov) + 'X_Coverage'].values:
                covMin = pol[str(min_cov) + 'X_Coverage'].value_counts()['Yes']/len(pol)*100
            else:
                covMin = 0
            if 'Yes' in pol['1000X_Coverage'].values:
                cov1000 = pol['1000X_Coverage'].value_counts()['Yes']/len(pol)*100
            else:
                cov1000 = 0

            # Calculate position diversity
            major = pol_nucl.max(axis=1)
            pol.insert(15, 'freq_major_base', major)
            pol.insert(16, 'theta', np.where((1-pol['freq_major_base']) > tc, 1, 0))
            pol.loc[pol['freq_major_base'] == 0, 'theta'] = 0
            nucl_diversity = pol_nucl * (1 - pol_nucl)
            pol.insert(17, 'diversity', nucl_diversity.sum(axis=1))

            # Calculate sample distance and Estimated time of infection (ETI)
            def calculate_ETI(pol, SID, covMin, cov1000, avg_cov, proj_name, eti_m, eti_c):

                #Calculate per base distance
                pol.insert(18, 'distance', pol['theta'] * pol['diversity'])

                # Filter out positions with too low base coverage
                pol_filtered = pol[pol[str(min_cov) + 'X_Coverage'] == 'Yes']

                # Calculate pairwise distance and ETI for sample
                pw_distance = pol_filtered['distance'].mean()
                pol.insert(19, 'avg_pairwise_distance', pw_distance)
                eti = eti_m * pw_distance + eti_c
                pol.insert(20, 'ETI', eti)
                
                # Print and return results
                sample_info = SID[:-25].split('_')
                #sample_info = sample_info[1:]   ##if sample name is prepended with 'X_'
                print(proj_name, sample_info[2], sample_info[0], pw_distance, eti, str(round(covMin, 2)) + '%', str(round(cov1000, 2)) + '%',  str(int(round(avg_cov))), sep=',')
                return pol

            pol = calculate_ETI(pol, SID, covMin, cov1000, avg_cov, proj_name, eti_m, eti_c)

            # Output
            pol.to_csv(outfile, sep=',', index=False)
