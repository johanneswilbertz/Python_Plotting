# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Generate dataframe
df = pd.read_csv(r'L:\ASD\Experiments\ASD_Exp11_D6_Pilot_Test_Dec_2020\Output_Data\CSV\processed_data_ASD_Exp11_D6_Pilot_Test_Dec_2020.csv') #_normed_by_dsmo.csv
# Drop column with just NANs or other undesired columns
df = df.drop(['Nuclei_dead'], axis=1)

# # Drop outliers based on z-score threshold
# def drop_numerical_outliers(df, z_thresh=4):
#     # Constrains will contain `True` or `False` depending on if it is a value below the threshold.
#     constrains = df.select_dtypes(include=[np.number]) \
#         .apply(lambda x: np.abs(stats.zscore(x, nan_policy='omit')) < z_thresh) \
#         .all(axis=1)
#     # Drop (inplace) values set to be rejected
#     df.drop(df.index[~constrains], inplace=True)
    
# drop_numerical_outliers(df)

# # Generation of dataframes to plot all HTs or WTs togehter, or just a single compound concentration

# # Select only rows fullfiling certain criteria

df['genotype_CPD'] = df['Label'] + "\n" + df['CPD_ID']

ht_dmso = df[(df['Experiment'].str.contains('HT')) & (df['Experiment'].str.contains('DMSO'))]
wt_dmso = df[(df['Experiment'].str.contains('WT')) & (df['Experiment'].str.contains('DMSO'))]
wt_ht_dmso = wt_dmso.append(ht_dmso)

wt_ht_dmso_licl = df[((df['Experiment'].str.contains('WT')) | (df['Experiment'].str.contains('HT'))) & ((df['Experiment'].str.contains('DMSO')) | (df['Experiment'].str.contains('LiCl')))]
wt_ht_dmso_ha1077 = df[((df['Experiment'].str.contains('WT')) | (df['Experiment'].str.contains('HT'))) & ((df['Experiment'].str.contains('DMSO')) | (df['Experiment'].str.contains('HA1077')))]

# ht_licl = df[(df['Experiment'].str.contains('HT')) & (df['Experiment'].str.contains('LiCL'))]
# wt_licl = df[(df['Experiment'].str.contains('WT')) & (df['Experiment'].str.contains('LiCL'))]

# # Append and adjust dataframe for plotting 
# ht = ht_dmso.append(ht_licl)
# ht['genotype'] = 'HT'
# wt = wt_dmso.append(wt_licl)
# wt['genotype'] = 'WT'
# wt_ht = ht.append(wt)
# # Creating extra column for global genotype combined with treatment condition
# wt_ht['genotype_CPD'] = wt_ht['genotype'] + "\n" + wt_ht['CPD_Concentration']
# # 2 mM LiCl only
# wt_ht_2mM = wt_ht[(wt_ht['Experiment'].str.contains('2 mM')) | (wt_ht['Experiment'].str.contains('DMSO'))]

# What to plot
data = wt_ht_dmso
xgrouping = 'genotype_CPD'
datacolumn = "avg_intensity_hucd"
ytitle = 'HuC/D intensity'
xtitle = 'HA1077 (uM)'

# # ------------ PLOTTING OF BARPLOT OVERLAYED WITH SWARMPLOT -------------

# # Overlay swarm and boxplot for CATEGORICAL DATA WITHOUT MULTIPLE PERTUBATIONS (for example WT/MUT comparison and single treatment)

# plt.figure()
# fig, ax = plt.subplots(figsize=(4, 6))
# ax = sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 3})
# ax = sns.boxplot(x=xgrouping, y=datacolumn, data=data, width=0.4, showfliers=False) #width=0.6 (optional)
# # Transparancy of boxplot filling color
# for patch in ax.artists:
#     r, g, b, a = patch.get_facecolor()
#     patch.set_facecolor((r, g, b, .2))
# ax = sns.stripplot(x=xgrouping, y=datacolumn, data=data, size=10, alpha=0.3)
# ax.set(xlabel=None)
# ax.set(ylabel=ytitle)

# ------------ PLOTTING OF LINEPLOT -------------

# # Line plot with error bars for CATEGORICAL DATA WITH MULTIPLE PERTUBATIONS (for example WT/MUT treated with multiple compound concentrations)

# plt.figure()
# fig, ax = plt.subplots(figsize=(10, 6))
# ax = sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 4, 'lines.markersize': 15})
# ax = sns.lineplot(data=data, x="Concentration", y=datacolumn, hue="Label", marker="o", err_style="bars", ci=99.7)
# ax.set(xlabel=xtitle, ylabel=ytitle)
# # Put the legend out of the figure and remove title
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles=handles[1:], labels=labels[1:], bbox_to_anchor=(1.02, 0.3), loc=3, borderaxespad=0., frameon=False)
# plt.tight_layout()

# # ------------ PLOTTING OF FACETGRIDPLOT (some elements not necessary with Seaborn catplot) -------------

# plt.figure()

# # Plot design
# sns.set_context("talk")
# sns.set_style("ticks")

# # Initialize a grid of plots with an Axes for each Label (not necessary for sns.catplot)
# # grid = sns.FacetGrid(df, col="Label", hue="Label", palette="pastel",
# #                       col_wrap=3, height=3)

# # Draw a line plot to show all the labels (col_wrap determines organization)
# grid = sns.catplot("Concentration", datacolumn, col='Label', hue="Label", data=df, kind="point", col_wrap=3, height=3, ci=68)

# # Draw a horizontal line to show the starting point
# grid.map(plt.axhline, y=1, ls=":", c=".5")

# # Adjust the tick positions and labels (not necessary for sns.catplot)
# # grid.set(xticks=np.arange(3), yticks=[0.5, 1, 1.5, 2],
# #           xlim=(-.1, 2.1), ylim=(0.4, 2.1))

# #Adjust axis labels
# grid.set_axis_labels(xtitle, ytitle)

# ------------ PLOTTING OF MULTIPLE PLOTS WITH DIFFERENT Y-AXIS -------------

plt.figure()
 
df1 = data.select_dtypes([np.number])
df1 = df1.drop(['Concentration'], axis=1)

df1_col = df1.columns

sns.set_context("talk")
sns.set_style("ticks")   

fig, axes = plt.subplots(nrows=5, ncols=7, figsize = (20, 15)) 

for idx, (col, ax) in enumerate(zip(df1.columns, axes.flatten())):
#for ax, feature in zip(axes.flatten(), df1.columns):
    sns.boxplot(ax=ax, x=xgrouping, y=df1[col], data=data, width=0.5, showfliers=False) #width=0.6 (optional)
    # Transparancy of boxplot filling color
    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .2))
    sns.stripplot(ax=ax, x=xgrouping, y=df1[col], data=data, size=10, alpha=0.3)
    plt.subplots_adjust(wspace=.5, hspace=.5) 
    ax.set(xlabel=None)  
else:
    [ax.set_visible(False) for ax in axes.flatten()[idx+1:]]

plt.tight_layout()
# plt.show()

# # ------------ PLOTTING OF MULTIPLE PLOTS WITH DIFFERENT Y-AXIS (LINEPLOT) -------------

# plt.figure()

# df1 = data.select_dtypes([np.int, np.float])
# df1 = df1.drop(['Concentration'], axis=1)

# sns.set_context("talk")
# sns.set_style("ticks")   

# fig, axes = plt.subplots(nrows=7, ncols=5, figsize = (30, 30)) 

# for idx, (col, ax) in enumerate(zip(df1.columns, axes.flatten())):
# #for ax, feature in zip(axes.flatten(), df1.columns):
#     sns.lineplot(ax = ax, data=data, x="Concentration", y=df1[col], hue='Label',marker="o", err_style="bars", ci=99.7)
#     plt.subplots_adjust(wspace=.5, hspace=.5) 
#     ax.get_legend().remove()
#     handles, labels = ax.get_legend_handles_labels()
#     ax.set(xlabel=xtitle)
    
# else:
#     [ax.set_visible(False) for ax in axes.flatten()[idx+1:]]

# ax.legend(handles=handles[1:], labels=labels[1:], bbox_to_anchor=(1.02, 0.3), loc=3, borderaxespad=0., frameon=False)
    
# plt.tight_layout()
# # plt.show()