# -*- coding: utf-8 -*-
"""
Johannes Wilbertz, PhD (Ksilink)

- Script to plot catergorical data as boxplots including all datapoints as swarmplot.
- Outliers can be removed based on Z-score.(based on: https://stackoverflow.com/a/56725366)
- Statistical comparisons can be added between groups. Different tests available (based on Statannot: https://github.com/webermarcolivier/statannot) 

"""

import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statannot import add_stat_annotation

# Generate dataframe
df = pd.read_csv(r'L:\ASD\Checkout_Results\210708_SynapsesD30_Drugs\agASD-60X-260521-Exp26-Hit-MAPK-pl1-Synapsin-MAP2-Homer-AT.csv') #_normed

# Data wrangling and renaming
df.loc[df['tags'].str.contains('HT3;LiCl'), 'tags'] = 'LiCl;HT3'
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

# What to plot
data = df
xgrouping = 'tags'
hue = "Plate"
datacolumn = 'synapses_per_nuclei'
medianline = 'DMSO;HT3'
ytitle = datacolumn
#xtitle = 'LiCl (mM)'
order=['DMSO;WT', 'LiCl;WT', 'GSK;WT', 'C401;WT', 'B-Raf1;WT', 'DMSO;HT3', 'LiCl;HT3', 'GSK;HT3', 'C401;HT3', 'B-Raf1;HT3']
stat_pairs = [("DMSO;WT", "DMSO;HT3"), ('DMSO;WT', "LiCl;WT"), ('DMSO;WT', "GSK;WT"), ('DMSO;WT', "C401;WT"), ('DMSO;WT', "B-Raf1;WT"),
              ('DMSO;HT3', "LiCl;HT3"), ('DMSO;HT3', "GSK;HT3"), ('DMSO;HT3', "C401;HT3"), ('DMSO;HT3', "B-Raf1;HT3")]
stat_test = 'Mann-Whitney'  # Test value should be one of the following: t-test_ind, t-test_welch, t-test_paired, 
                            # Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal
                            
                            # Mann-Whitney can be used to compare differences between two independent groups when the dependent variable 
                            # is either ordinal or continuous, but not normally distributed. Good non-parametric alternative for t-test.
                            

# ------------ PLOTTING OF BARPLOT OVERLAYED WITH SWARMPLOT -------------

# Overlay swarm and boxplot for CATEGORICAL DATA WITHOUT MULTIPLE PERTUBATIONS (for example WT/MUT comparison and single treatment)

ax = sns.set_context("talk")
ax = sns.boxplot(x=xgrouping, y=datacolumn, data=data, width=0.7, showfliers=False, zorder=1,
                 order=order)
# Transparancy of boxplot filling color
for patch in ax.artists:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, 0))
ax = sns.stripplot(x=xgrouping, y=datacolumn, data=data, size=10, alpha=0.3, zorder=0, hue=xgrouping,
                   order=order)
#Figure aesthetics
ax.set(xlabel=None)
ax.set(ylabel=ytitle)
sns.despine()
ax.legend_.remove() 
# ax.set(ylim=(-0.2, 60))
ax.set_xticklabels(ax.get_xticklabels(),rotation=45, ha='right', rotation_mode="anchor")
# # ax.set(yscale="log")

# Calculation of mean dashed line for easier comparison
medianline_data = df[(df[xgrouping].str.contains(medianline))]
median = medianline_data[datacolumn].median()
ax.axhline(y=median, ls=":", c=".5")

# ------------------------------------------------------------------------

# Add statistical annotation on top of figure

# OPTION 1: Stats plotting for total population of all technical replicates from all biological replicates. Careful: High N = small p-value. Look at effect size!
add_stat_annotation(ax, data=data, x=xgrouping, y=datacolumn, order=order,
                    box_pairs=stat_pairs, test=stat_test, text_format='star', loc='outside', verbose=2)

# # OPTION 2: Stats plotting only for means of biological replicates.
# ReplicateAverages = data.groupby([xgrouping, hue], as_index=False).agg({datacolumn: "median"})
# ReplicateAverages = ReplicateAverages.sort_values(['tags', 'Plate']) # Sorting is necessary to plot averages accurately on dataclouds
# add_stat_annotation(ax, data=ReplicateAverages, x=xgrouping, y=datacolumn, order=["GIBCO control", "C89 control", "SNCA triplication"],
#                     box_pairs=[("GIBCO control", "C89 control"), ("GIBCO control", "SNCA triplication"), ("C89 control", "SNCA triplication")],
#                     test='t-test_welch', text_format='star', loc='outside', verbose=2)



