# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 16:15:05 2022

@author: JohannesWilbertz
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

df  = pd.read_csv('L:\\PROJECTS\\PD\\Experiment\\PrimaryScreen\\CSV\\df_lda3_classes_all_normed_data.csv')

# Dropping the rows of unwanted plates and fluorescent compounds
df = df[df['Plate'].str.contains('PD-SC2-01') == False]
df = df[df['tags'].str.contains('Terthiophene|Rhodamine B|Obatoclax Mesylate') == False]

# Group by tag and average values
#df_grouped = df.groupby(['tags'], as_index=False).median()
#df_grouped.to_csv('L:\\PROJECTS\\PD\\Experiment\\PrimaryScreen\\CSV\\df_lda3_classes_all_normed_data_GROUPED.csv')

# Extract DMSO controls only
df_mut_dmso = df[df['tags'].str.contains('Mut;DMSO')]

# Calculate DMSO medians per plate 
df_median = df_mut_dmso.groupby(df_mut_dmso['Plate']).median().reset_index().add_suffix('_median')
df_median = df_median.rename(columns={'Plate_median': "Plate"})

# Calculate DMSO SDs per plate 
df_std = df_mut_dmso.groupby(df_mut_dmso['Plate']).std().reset_index().add_suffix('_std')
df_std = df_std.rename(columns={'Plate_std': "Plate"})

# Merge median and std values
df_median_std = pd.merge(df_median[['Plate','Nuclei_Number_Living_median', 'Cell_Neurites_LengthPerNuclei_MAP2_median']],
                         df_std[['Plate','Nuclei_Number_Living_std', 'Cell_Neurites_LengthPerNuclei_MAP2_std']],
                         on='Plate', how='left')
# Join medians and SDs per plate with all data
df = df.join(df_median_std.set_index('Plate'), on='Plate')

# Thresholds for toxcity selection
SD_nuclei_thresh = 3
SD_neurites_thresh = 3 

# Calculate toxicity limits
limit_nuclei = df['Nuclei_Number_Living_median']-SD_nuclei_thresh*df['Nuclei_Number_Living_std'] 
limit_neurites = df['Cell_Neurites_LengthPerNuclei_MAP2_median']-SD_neurites_thresh*df['Cell_Neurites_LengthPerNuclei_MAP2_std']

# Apply toxicity limits to dataset
df_non_toxic = (df['Nuclei_Number_Living']>=limit_nuclei) & (df['Cell_Neurites_LengthPerNuclei_MAP2']>=limit_neurites)
df = df.assign(non_toxic=df_non_toxic)

# Add label
df['exp_condition'] = np.where(df['tags'].isin(['Mut;DMSO','WT;DMSO','Mut;PRO']), df['tags'], 'Mut;Compound')
  
# Save as .csv
df.to_csv('L:\\PROJECTS\\PD\\Experiment\\PrimaryScreen\\CSV\\df_lda3_classes_all_normed_data_NUCLEI_NEURITE_MEDIANS.csv')

# Select only non toxic molecules
df_compounds = df[df.non_toxic]

# HIT SELECTION
# Remove controls
Search_for_These_values = ['DMSO','PRO'] 
pattern = '|'.join(Search_for_These_values)
df_compounds_noCtrls = df_compounds.loc[~df_compounds['tags'].str.contains(pattern, case=False)]
# Choose hits based on threshold
hit_selection_thresh = 3.5
df_pre_hits = df_compounds_noCtrls.loc[df_compounds_noCtrls['mahalanobis_from_neg'] > hit_selection_thresh]

# Drop all tags that occur only once
df_hits = df_pre_hits[df_pre_hits.duplicated(subset=['tags'], keep=False)]
df_hits_grouped = df_hits.groupby(['tags'], as_index=False).mean()

# Save individual replicates and mean of replicates as .csv
df_hits.to_csv('L:\\PROJECTS\\PD\\Experiment\\PrimaryScreen\\CSV\\df_lda3_classes_all_normed_data_ALL_HITS.csv')
df_hits_grouped.to_csv('L:\\PROJECTS\\PD\\Experiment\\PrimaryScreen\\CSV\\df_lda3_classes_all_normed_data_MEAN_HITS.csv')

# PLOT DATA
# Define compound of interest
compound_of_interest = 'Mut;WAY-657095'
df_compounds['exp_condition'] = np.where((df_compounds.tags == compound_of_interest),compound_of_interest,df_compounds.exp_condition)

sns.set_context("talk")
sns.set_style("ticks")

# Set colors
color_dict = dict({'WT;DMSO':'blue',
                  'Mut;Compound':'gray',
                  'Mut;DMSO':'red',
                  'Mut;PRO': 'orange',
                  compound_of_interest: 'black'})

hue_order=['WT;DMSO', 'Mut;PRO', 'Mut;Compound', 'Mut;DMSO', compound_of_interest]
ax = sns.scatterplot(data=df_compounds.sort_values('exp_condition', key=np.vectorize(hue_order.index)), 
                     x='mahalanobis_from_neg' ,y='mahalanobis_from_pos', 
                     hue='exp_condition', s=30,
                     edgecolor='none', hue_order=hue_order, palette=color_dict)

# Figure aesthetics 
sns.despine()
#ax.set_xlim(-0.01, 10)
#ax.set_ylim(-0.01, 30)
#ax.set(xlabel='Ratio HuC/D+ nuclei')
#ax.set(ylabel='Ratio Ki67+ nuclei')
#ax.axhline(y=limit_ki67, ls=":", c="0.5")
ax.axvline(x=hit_selection_thresh, ls=":", c="0.5")
plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0., frameon=False)