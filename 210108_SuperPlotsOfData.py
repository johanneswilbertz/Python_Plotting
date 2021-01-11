# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 17:33:17 2021

@author: Johannes Wilbertz, Ksilink

This script apends and visualizes the results of multiple imaging data files, representing biological replicates.
Also the technical replicates per biological replicate are plotted and p-value statistics are calculated based on the 
biological replicates. 

Based on idea & Python code provided here: https://rupress.org/jcb/article/219/6/e202001064/151717/SuperPlots-Communicating-reproducibility-and

"""
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Generate dataframes containing individual replicates
df1 = pd.read_csv(r'L:\PD\Experiment\Maturation_18\Output_Data\CSV\raw_data_normed_median_dmso_DRC_Tool.csv')
df2 = pd.read_csv(r'L:\PD\Experiment\Maturation_20\Output_Data\CSV\raw_data_normed_median_dmso_with_conditions.csv')

# Add plate/replicate number
df1.loc[df1['Plate'].str.contains('p1'), 'Plate'] = 'Plate1'
df1.loc[df1['Plate'].str.contains('p2'), 'Plate'] = 'Plate2'
df2.loc[df2['Plate'].str.contains('p1'), 'Plate'] = 'Plate3'
df2.loc[df2['Plate'].str.contains('p2'), 'Plate'] = 'Plate4'

# Append dataframes and adjust organisation for plotting
combined = df1.append(df2)
combined.loc[combined['CPD_ID'].str.contains('Prostratine'), 'CPD_ID'] = 'Prostratin' # Correct misspelling
combined['Class_CPD'] = combined['Class'] + ' ' + combined['CPD_ID']
combined = combined.sort_values(['Class_CPD', 'Plate']) # Sorting is necessary to plot averages accurately on dataclouds
combined_columns = combined.columns

# Drop outliers based on z-score threshold
def drop_numerical_outliers(combined, z_thresh=3):
    # Constrains will contain `True` or `False` depending on if it is a value below the threshold.
    constrains = combined.select_dtypes(include=[np.number]) \
        .apply(lambda x: np.abs(stats.zscore(x, nan_policy='omit')) < z_thresh) \
        .all(axis=1)
    # Drop (inplace) values set to be rejected
    combined.drop(combined.index[~constrains], inplace=True)
    
drop_numerical_outliers(combined)

# Data for plotting and statistics calculations
xgrouping = "Class_CPD" # Treatment Category
xgroup1 = "WT DMSO" # Statistics compared to this group
replicate = "Plate" #Replicate
datacolumn = "avg_intensity_TH_SNCA_SNCA"
ytitle = 'Norm.' + ' ' + datacolumn
plot_order=["WT DMSO", "WT PEP005", "WT Prostratin", "Mut DMSO", "Mut PEP005", "Mut Prostratin"]

# Generating a SuperPlotOfData
plt.figure()
fig, ax = plt.subplots(figsize=(10, 6))
ax = sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 3})

ReplicateAverages = combined.groupby([xgrouping, replicate], as_index=False).agg({datacolumn: "median"})
ReplicateAverages = ReplicateAverages.sort_values(['Class_CPD', 'Plate']) # Sorting is necessary to plot averages accurately on dataclouds
ReplicateAvePivot = ReplicateAverages.pivot_table(columns=xgrouping, values=datacolumn, index=replicate) 

ax = sns.stripplot(x=xgrouping, y=datacolumn, hue=replicate, size=7, data=combined, order=plot_order) 
ax = sns.stripplot(x=xgrouping, y=datacolumn, hue=replicate, size=20, edgecolor="k", linewidth=2, data=ReplicateAverages, order=plot_order)

# Figure aesthetics
ax.grid(False)
ax.legend_.remove() 
ax.set(xlabel=None)
ax.set(ylabel=ytitle)
ax.set_ylim([0, combined[datacolumn].max()])
ax.set_xticklabels(ax.get_xticklabels(),rotation=30)
sns.despine()

# Stats plotting based on ordered groups
for g in range(1, len(plot_order)): 
    xgroup2 = (plot_order[g]) 
    statistic, pvalue = scipy.stats.ttest_rel(ReplicateAvePivot[xgroup1], ReplicateAvePivot[xgroup2])
    P_value = str(float(round(pvalue, 3))) 
    x1, x2 = 0, g
    y, h, col = combined[datacolumn].max() + (g-1)/14, (g-1)/14, 'k'
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col) 
    ax.text((x1+x2)*.5, y+h, "P = "+P_value, ha='center', va='bottom', color=col, size=18)
