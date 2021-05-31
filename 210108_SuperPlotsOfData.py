# -*- coding: utf-8 -*-
"""
Creating SuperPlotsOfData with statistics --> Overlay if two swarmplots creates transparent
display of all underlying data, including technical and biological replicates.

Johannes Wilbertz, PhD @ Ksilink, 31 May 2020

Credits:

1. SuperPlotsOfData is a concept to visualize technical and biological replicates transpaently 
and was publsihed by Lord et al. JCB, 2020 (https://doi.org/10.1083/jcb.202001064). 

2. The package "statanot" for plotting statistics outside of the figure is form webermarcolivier 
and colleagues (https://github.com/webermarcolivier/statannot).

3. The code for removal of outliers in numerical and non-numerical datasets is originally from
KeyMaker00 (https://stackoverflow.com/a/56725366).

"""
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statannot import add_stat_annotation

# Generate dataframes containing individual replicates
df = pd.read_csv(r'L:\PD\Experiment\Maturation_36\Output_Data\CSV\raw_data_normed_median_dmso_snca.csv')

# Data wrangling if needed

# # Drop columns which are not needed
# df = df.drop(columns=['Nuclei_Big'])

# # Rename some data tags for consistency
# df.loc[df['tags'].str.contains('TEST123'), 'tags'] = 'TESTABC' 
# df.loc[df['tags'].str.contains('TEST456'), 'tags'] = 'TESTDEF' 

# Adjust organisation for plotting
df = df.sort_values(['tags', 'Plate']) # Sorting is necessary to plot averages accurately on dataclouds
df_columns = df.columns
combined = df

# # Drop outliers based on z-score threshold
# def drop_numerical_outliers(combined, z_thresh=3):
#     # Constrains will contain `True` or `False` depending on if it is a value below the threshold.
#     constrains = combined.select_dtypes(include=[np.number]) \
#         .apply(lambda x: np.abs(stats.zscore(x, nan_policy='omit')) < z_thresh) \
#         .all(axis=1)
#     # Drop (inplace) values set to be rejected
#     combined.drop(combined.index[~constrains], inplace=True)
    
# drop_numerical_outliers(combined)

# Data for plotting and statistics calculations
xgrouping = "tags" # Treatment Category
replicate = "Plate" #Kind of used replicate
datacolumn = "avg_intensity_SNCA" # Nuclei_Tot, Nuclei_dead, avg_intensity_SNCA, sum_intensity_per_nuclei_SNCA
                                   # ratio_nuclei_MAP2, surface_MAP2, branching_points_MAP2_per_nuclei, skelet_length_MAP2_per_nuclei
ytitle = 'Norm.' + ' ' + datacolumn
plot_order=["Ctrl;DMSO", "Ctrl;PEP005", "Ctrl;Prostratin", "Tripli;DMSO", "Tripli;PEP005", "Tripli;Prostratin"]

# Generating a SuperPlotOfData

# Create new figure and two subplots
ax = sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 3})

ReplicateAverages = combined.groupby([xgrouping, replicate], as_index=False).agg({datacolumn: "median"})
ReplicateAverages = ReplicateAverages.sort_values(['tags', 'Plate']) # Sorting is necessary to plot averages accurately on dataclouds

ax = sns.stripplot(x=xgrouping, y=datacolumn, hue=replicate, size=7, data=combined, order=plot_order) 
ax = sns.stripplot(x=xgrouping, y=datacolumn, hue=replicate, size=20, edgecolor="k", linewidth=2, data=ReplicateAverages, order=plot_order)

# Figure aesthetics
ax.grid(False)
ax.legend_.remove() 
ax.set(xlabel=None)
ax.set(ylabel=ytitle)
ax.set_xticklabels(ax.get_xticklabels(),rotation=45, ha='right', rotation_mode="anchor")
sns.despine()

# # OPTION 1: Stats plotting for total population of all technical replicates from all biological replicates. Careful: High N = small p-value. Look at effect size!
# add_stat_annotation(ax, data=combined, x=xgrouping, y=datacolumn, order=plot_order,
#                     box_pairs=[("Ctrl;DMSO", "Tripli;DMSO"), ("Ctrl;DMSO", "Ctrl;PEP005"), ("Tripli;DMSO", "Tripli;PEP005")],
#                     test='Mann-Whitney', text_format='star', loc='outside', verbose=2)

# OPTION 2: Stats plotting only for means of biological replicates.   
add_stat_annotation(ax, data=ReplicateAverages, x=xgrouping, y=datacolumn, order=plot_order,
                    box_pairs=[("Ctrl;DMSO", "Tripli;DMSO"), ("Ctrl;DMSO", "Ctrl;PEP005"), ("Tripli;DMSO", "Tripli;PEP005")],
                    test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
