# -*- coding: utf-8 -*-
"""
Johannes Wilbertz, PhD (Ksilink)

- Script to plot catergorical data as boxplots including all datapoints as swarmplot.
- Outliers can be removed based on Z-score.(based on: https://stackoverflow.com/a/56725366)
- Statistical comparisons can be added between groups. Different tests available (based on Statannot: https://github.com/webermarcolivier/statannot) 

"""

import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from statannot import add_stat_annotation

# Plot all image features ('yes') or only select ones ('no')?
plot_all_features = 'no'

# Specify input & generate dataframe
df = pd.read_csv(r'L:\PD\Experiment\Maturation_49\ICC\CSV\agSNCA-Mat49-x40-ICC-MAP2-SNCA-TH-AV.csv')
if plot_all_features == 'yes':
    results_dir = 'L:\\PD\Experiment\\Maturation_49\\ICC\Graphs\\' # use double backslash, otherwise interpreted as "escape"
    plt.rcParams["savefig.directory"] = os.chdir(os.path.dirname(results_dir))
df = df[df['tags'].notna()]
df2 = df.copy()

# Get unique data tags/categories
unique_tags = df['tags'].unique()
print(unique_tags)

# # Data wrangling and renaming if needed
# df['tags'] = df['Condition'] + ' ' + df['Concentration'] + ' ' + df['CPD_ID']
# df.loc[df['tags'].str.contains('WT2 0µM DMSO'), 'tags'] = 'WT2 DMSO'
# df.loc[df['tags'].str.contains('HT3 0µM DMSO'), 'tags'] = 'HT3 DMSO'
# df2 = df.copy()

# Filter outliers based on sigma threshold per data category (IF DIFFERENT PLATES ARE USED  NORMALIZE THEM FIRST)
list_df_tag = []
for tag in unique_tags:
    df_tag = df.loc[df['tags'] == tag]
    
    # Dataframe is pre-processed to exclude inf values or columns with only single value since this won't allow z-score calculation
    df_tag.replace([np.inf, -np.inf], np.nan, inplace=True)     # Change all inf for NaN
    df_numbers = df_tag.select_dtypes(include=[np.number])      # Choose only numerical values
    nunique = df_numbers.nunique()                              # Number of unique values per column
    cols_to_drop = nunique[nunique < 2].index                   # Identify columns with only 1 value (= same number) or 0 values (= all NaN)
    df_tag = df_tag.drop(cols_to_drop, axis=1)                  # Drop all these columns from main dataframe 
    
    # Drop outliers based on z-score threshold
    def drop_numerical_outliers(df_tag, z_thresh=3):
        # Constrains will contain `True` or `False` depending on if it is a value below the threshold.
        constrains = df_tag.select_dtypes(include=[np.number]) \
            .apply(lambda x: np.abs(stats.zscore(x, nan_policy='omit')) < z_thresh) \
            .all(axis=1)
        # Drop (inplace) values set to be rejected
        df_tag.drop(df_tag.index[~constrains], inplace=True)
        
    drop_numerical_outliers(df_tag)
    list_df_tag.append(df_tag)
    
df = pd.concat(list_df_tag)

# Which columns are left after outlier removal?
df_columns = df.columns
                         
# What to plot?
data = df
xgrouping = 'tags'
#hue = "Plate"
datacolumn = 'Nuclei_Number_Living' # 'Cell_TH_SNCA_Intensity_SumIntensityPerNuclei_SNCA', Cell_TH_SNCA_Intensity_MeanIntensity_SNCA
medianline = 'WT DMSO'
ytitle = datacolumn
#xtitle = 'LiCl (mM)'
order= ['WT DMSO', 'WT Prostratin', 'Tripli DMSO', 'Tripli Prostratin']  
stat_pairs = [('WT DMSO', 'WT Prostratin'), ('WT DMSO', 'Tripli DMSO'),
              ('Tripli DMSO', 'Tripli Prostratin')]
stat_test = 'Mann-Whitney'  # Test value should be one of the following: t-test_ind, t-test_welch, t-test_paired, 
                            # Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal
                            
                            # Mann-Whitney can be used to compare differences between two independent groups when the dependent variable 
                            # is either ordinal or continuous, but not normally distributed. Good non-parametric alternative for t-test.

# ------------ PLOTTING OF BARPLOT OVERLAYED WITH SWARMPLOT -------------

# Overlay swarm and boxplot for CATEGORICAL DATA WITHOUT MULTIPLE PERTUBATIONS (for example WT/MUT comparison and single treatment)

pal = sns.color_palette()   # Get all color codes
pal = pal.as_hex()          # Transform color codes to hexformat
#sns.color_palette(['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd'])
# use colors in this order
palette = {'WT DMSO': '#1f77b4', 'WT Prostratin': '#ff7f0e', 
           'Tripli DMSO': '#2ca02c', 'Tripli Prostratin': '#d62728'} 

if plot_all_features == 'no':
    ax = sns.set_context("talk")
    ax = sns.boxplot(x=xgrouping, y=datacolumn, data=data, width=0.7, 
                     showfliers=False, order=order, zorder=1)
    # Transparancy of boxplot filling color
    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, 0))
    ax = sns.stripplot(x=xgrouping, y=datacolumn, data=data, size=10, alpha=0.3, 
                       zorder=0, hue=xgrouping, order=order, palette=palette)
    #Figure aesthetics
    ax.set(xlabel=None)
    ax.set(ylabel=ytitle)
    sns.despine()
    ax.legend_.remove() 
    # ax.set(ylim=(-0.2, 60))
    ax.set_xticklabels(ax.get_xticklabels(),rotation=45, ha='right', 
                       rotation_mode="anchor")
    # # ax.set(yscale="log")
    
    # Calculation of mean dashed line for easier comparison
    medianline_data = data[(data[xgrouping].str.contains(medianline))]
    median = medianline_data[datacolumn].median()
    ax.axhline(y=median, ls=":", c=".5")
    
    # Add statistical annotation on top of figure
    # OPTION 1: Stats plotting for total population of all technical replicates from all biological replicates. 
    # Careful: High N = small p-value. Look at effect size!
    add_stat_annotation(ax, data=data, x=xgrouping, y=datacolumn, order=order,
                    box_pairs=stat_pairs, test=stat_test, text_format='star', loc='outside', verbose=2)

if plot_all_features == 'yes':
    df_numeric = df.select_dtypes([np.int, np.float])
    for i, col in enumerate(df_numeric.columns):
        plt.figure(i)
        ax = sns.set_context("talk")
        ax = sns.boxplot(x=xgrouping, y=col, data=data, width=0.7, showfliers=False, zorder=1,
                         order=order)
        # Transparancy of boxplot filling color
        for patch in ax.artists:
            r, g, b, a = patch.get_facecolor()
            patch.set_facecolor((r, g, b, 0))
        ax = sns.stripplot(x=xgrouping, y=col, data=data, size=10, alpha=0.3, zorder=0, hue=xgrouping,
                       order=order, palette=palette)
        #Figure aesthetics
        ax.set(xlabel=None)
        ax.set(ylabel=col)
        sns.despine()
        ax.legend_.remove() 
        # ax.set(ylim=(-0.2, 60))
        ax.set_xticklabels(ax.get_xticklabels(),rotation=45, ha='right', rotation_mode="anchor")
        # # ax.set(yscale="log")
    
        # Calculation of mean dashed line for easier comparison
        medianline_data = data[(data[xgrouping].str.contains(medianline))]
        median = medianline_data[col].median()
        ax.axhline(y=median, ls=":", c=".5")
        
        # Add statistical annotation on top of figure
        # OPTION 1: Stats plotting for total population of all technical replicates from all biological replicates. 
        # Careful: High N = small p-value. Look at effect size!
        add_stat_annotation(ax, data=data, x=xgrouping, y=col, order=order,
                        box_pairs=stat_pairs, test=stat_test, text_format='star', loc='outside', verbose=2)
        
        plt.savefig(col + '.png', bbox_inches='tight')
# ------------------------------------------------------------------------

# Add statistical annotation on top of figure

# # OPTION 1: Stats plotting for total population of all technical replicates from all biological replicates. Careful: High N = small p-value. Look at effect size!
# add_stat_annotation(ax, data=data, x=xgrouping, y=datacolumn, order=order,
#                     box_pairs=stat_pairs, test=stat_test, text_format='star', loc='outside', verbose=2)

# # # OPTION 2: Stats plotting only for means of biological replicates.
# ReplicateAverages = data.groupby([xgrouping, hue], as_index=False).agg({datacolumn: "median"})
# ReplicateAverages = ReplicateAverages.sort_values(['tags', 'Plate']) # Sorting is necessary to plot averages accurately on dataclouds
# add_stat_annotation(ax, data=ReplicateAverages, x=xgrouping, y=datacolumn, order=["GIBCO control", "C89 control", "SNCA triplication"],
#                     box_pairs=[("GIBCO control", "C89 control"), ("GIBCO control", "SNCA triplication"), ("C89 control", "SNCA triplication")],
#                     test='t-test_welch', text_format='star', loc='outside', verbose=2)



