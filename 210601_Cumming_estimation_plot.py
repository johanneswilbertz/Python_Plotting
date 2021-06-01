# -*- coding: utf-8 -*-
"""
Johannes Wilbertz @ Ksilink, 01 June 2021

Based on:

Moving beyond P values: Everyday data analysis with estimation plots
Joses Ho, Tayfun Tumkaya, Sameer Aryal, Hyungwon Choi, Adam Claridge-Chang
Nature Methods 2019, 1548-7105. 10.1038/s41592-019-0470-3

"""
import pandas as pd
import dabest
import matplotlib.pyplot as plt

# Generate dataframes containing individual replicates
df = pd.read_csv(r'L:\PD\Experiment\Maturation_35\Output_Data\CSV_400\raw_data_normed_median_dmso_snca_400.csv')

# Rename some data tags for consistency
df.loc[df['tags'].str.contains('Prostratin;WT'), 'tags'] = 'WT+Pro' 
df.loc[df['tags'].str.contains('PEP005;WT'), 'tags'] = 'WT+PEP' 
df.loc[df['tags'].str.contains('WT;DMSO'), 'tags'] = 'WT+DMSO' 
df.loc[df['tags'].str.contains('Mut;DMSO'), 'tags'] = 'GS+DMSO' 
df.loc[df['tags'].str.contains('Mut;PEP005'), 'tags'] = 'GS+PEP' 
df.loc[df['tags'].str.contains('Mut;Prostratin'), 'tags'] = 'GS+Pro' 

# Overview of data categories
df_columns = df.columns

# Load the above data into `dabest`.
df_dabest = dabest.load(data=df, x="tags", y="Nuclei_dead",
                          idx=("WT+DMSO", "WT+PEP", "WT+Pro", "GS+DMSO", "GS+PEP", "GS+Pro"))

# Produce a Cumming estimation plot. Choose "median" instead of "mean" for skewed distributions. 
plt.style.use("seaborn-ticks")

f = df_dabest.median_diff.plot(raw_marker_size=3); #contrast_ylim=(-0.15, 0.4)

rawswarm_axes = f.axes[0]
contrast_axes = f.axes[1]

# rawswarm_axes.set_xticklabels(rawswarm_axes.get_xticklabels(),rotation=45, ha='right', rotation_mode="anchor")
# contrast_axes.set_xticklabels(contrast_axes.get_xticklabels(),rotation=45, ha='right', rotation_mode="anchor")


"""
The tutorial can be found here:
    
    https://acclab.github.io/DABEST-python-docs/tutorial.html
    
"""