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

# Generate dataframe
df = pd.read_csv(r'L:\PD\Experiment\Maturation_26\Output_Data\CSV\raw_data.csv') #_normed

# Data wrangling and renaming
df = df[df.tags.str.contains("no treat")] # Add ~df.tags.str.contains("no treat") to exclude, remove to include
df = df[~df.tags.str.contains("Edi001A1/2")]
df.loc[df['tags'].str.contains('A18944'), 'tags'] = 'GIBCO control'
df.loc[df['tags'].str.contains('c89bmS4'), 'tags'] = 'C89 control'
df.loc[df['tags'].str.contains('Edi001A'), 'tags'] = 'SNCA triplication'

# Adjust organisation for plotting
df = df.sort_values(['tags', 'Plate']) # Sorting is necessary to plot averages accurately on dataclouds
df_columns = df.columns

# Drop outliers based on z-score threshold
# Dataframe is pre-processed to exclude inf values or columns with only single value since this won't allow z-score calculation
df.replace([np.inf, -np.inf], np.nan, inplace=True)     # Change all inf for NaN
df_numbers = df.select_dtypes(include=[np.number])      # Choose only numerical values
nunique = df_numbers.nunique()                          # Number of unique values per column
cols_to_drop = nunique[nunique < 2].index               # Identify columns with only 1 value (= same number) or 0 values (= all NaN)
df = df.drop(cols_to_drop, axis=1)                      # Drop all these columns from main dataframe
def drop_numerical_outliers(df, z_thresh=3):
    # Constrains will contain `True` or `False` depending on if it is a value below the threshold.
    constrains = df.select_dtypes(include=[np.number]) \
        .apply(lambda x: np.abs(stats.zscore(x, nan_policy='omit')) < z_thresh) \
        .all(axis=1)
    # Drop (inplace) values set to be rejected
    df.drop(df.index[~constrains], inplace=True)
    
drop_numerical_outliers(df)

# Data for plotting and statistics calculations
xgrouping = "tags" # Treatment Category
replicate = "Plate" #Kind of used replicate
datacolumn = "branching_points_MAP2_per_nuclei" # Nuclei_Tot, Nuclei_dead, avg_intensity_SNCA, sum_intensity_per_nuclei_SNCA
                                   # ratio_nuclei_MAP2, surface_MAP2, branching_points_MAP2_per_nuclei, skelet_length_MAP2_per_nuclei
ytitle = datacolumn
plot_order=["GIBCO control", "C89 control", "SNCA triplication"]
stat_pairs = [("GIBCO control", "C89 control"), ("GIBCO control", "SNCA triplication"), ("SNCA triplication", "C89 control")]
test='t-test_welch' # Test value should be one of the following: t-test_ind, t-test_welch, t-test_paired, 
                    # Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal
                    # Mann-Whitney can be used to compare differences between two independent groups when the dependent variable 
                    # is either ordinal or continuous, but not normally distributed. Good non-parametric alternative for t-test.     

# Generating a SuperPlotOfData

# Create new figure and two subplots
ax = sns.set_context("talk", font_scale=1.3) #rc={"lines.linewidth": 3}
ReplicateAverages = df.groupby([xgrouping, replicate], as_index=False).agg({datacolumn: "median"})
ReplicateAverages = ReplicateAverages.sort_values(['tags', 'Plate']) # Sorting is necessary to plot averages accurately on dataclouds

ax = sns.stripplot(x=xgrouping, y=datacolumn, hue=replicate, size=7, data=df, order=plot_order) 
ax = sns.stripplot(x=xgrouping, y=datacolumn, hue=replicate, size=20, edgecolor="k", linewidth=2, data=ReplicateAverages, order=plot_order)

# Figure aesthetics
ax.grid(False)
ax.legend_.remove() 
ax.set(xlabel=None)
ax.set(ylabel=ytitle)
ax.set_xticklabels(ax.get_xticklabels(),rotation=45, ha='right', rotation_mode="anchor")
sns.despine()

# OPTION 1: Stats plotting for total population of all technical replicates from all biological replicates. Careful: High N = small p-value. Look at effect size!
add_stat_annotation(ax, data=df, x=xgrouping, y=datacolumn, order=plot_order,
                    box_pairs=stat_pairs,
                    test=test, text_format='star', loc='outside', verbose=2)

# # OPTION 2: Stats plotting only for means of biological replicates.   
# add_stat_annotation(ax, data=ReplicateAverages, x=xgrouping, y=datacolumn, order=plot_order,
#                     box_pairs=stat_pairs,
#                     test=test, text_format='star', loc='outside', verbose=2)
