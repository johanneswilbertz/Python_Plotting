# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 11:13:10 2020

@author: JohannesWilbertz
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import math

# Generate dataframe
df = pd.read_csv('L:\ASD\Experiments\ASD_Exp8_D14_LiCl_Nov_2020\Output_Data\CSV\processed_data_ASD_Exp8_D14_LiCl_Nov_2020_normed_by_dsmo.csv') #_normed_by_dsmo.csv
# Drop column with just NANs or other undesired columns
df = df.drop(['Nuclei_dead'], axis=1)

# Drop outliers based on z-score threshold
def drop_numerical_outliers(df, z_thresh=4):
    # Constrains will contain `True` or `False` depending on if it is a value below the threshold.
    constrains = df.select_dtypes(include=[np.number]) \
        .apply(lambda x: np.abs(stats.zscore(x)) < z_thresh) \
        .all(axis=1)
    # Drop (inplace) values set to be rejected
    df.drop(df.index[~constrains], inplace=True)
    
drop_numerical_outliers(df)

# GroupBy, reset index and concatanate two dataframes
grouped_df = df.groupby(["Label", "Concentration"]).mean()
grouped_df = grouped_df.reset_index(level=['Label', 'Concentration'])
grouped_df_std = df.groupby(["Label", "Concentration"]).std()
grouped_df_std = grouped_df_std.reset_index(level=['Label', 'Concentration'])
grouped_df_std.columns = [str(col) + '_std' for col in grouped_df_std.columns]
grouped_df = pd.concat([grouped_df, grouped_df_std.reindex(grouped_df.index)], axis=1)

# Select only rows fullfiling certain criteria
ht_dmso = df[(df['Experiment'].str.contains('HT')) & (df['Experiment'].str.contains('DMSO'))]
ht_licl = df[(df['Experiment'].str.contains('HT')) & (df['Experiment'].str.contains('LiCL'))]

wt_dmso = df[(df['Experiment'].str.contains('WT')) & (df['Experiment'].str.contains('DMSO'))]
wt_licl = df[(df['Experiment'].str.contains('WT')) & (df['Experiment'].str.contains('LiCL'))]
# Append and adjust dataframe for plotting 
ht = ht_dmso.append(ht_licl)
ht['genotype'] = 'HT'
wt = wt_dmso.append(wt_licl)
wt['genotype'] = 'WT'
wt_ht = ht.append(wt)
# Creating extra column for global genotype combined with treatment condition
wt_ht['genotype_CPD'] = wt_ht['genotype'] + "\n" + wt_ht['CPD_Concentration']
# 2 mM LiCl only
wt_ht_2mM = wt_ht[(wt_ht['Experiment'].str.contains('2 mM')) | (wt_ht['Experiment'].str.contains('DMSO'))]

ht3 = df[(df['Experiment'].str.contains('HT3'))]

# ------------ PLOTTING MULTIPLE SUBPLOTS -------------

df1 = wt_ht.select_dtypes([np.int, np.float])

sns.set_context("talk")
sns.set_style("ticks")   

fig, axes = plt.subplots(nrows=7, ncols=5, figsize = (30, 30), sharex=True)

for idx, (col, ax) in enumerate(zip(df1.columns, axes.flatten())):
#for ax, feature in zip(axes.flatten(), df1.columns):
    sns.lineplot(ax = ax, data=wt_ht, x="Concentration", y=df1[col], hue='genotype',marker="o", err_style="bars", ci=68)
    plt.subplots_adjust(wspace=.5, hspace=.5) 
    ax.get_legend().remove()
    handles, labels = ax.get_legend_handles_labels()
    
else:
    [ax.set_visible(False) for ax in axes.flatten()[idx+1:]]

ax.legend(handles=handles[1:], labels=labels[1:], bbox_to_anchor=(1.02, 0.3), loc=3, borderaxespad=0., frameon=False)
#fig.legend(handles, labels, loc='lower right', title=None)
    

plt.tight_layout()
# plt.show()