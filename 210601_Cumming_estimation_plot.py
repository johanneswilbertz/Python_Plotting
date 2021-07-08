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
df = pd.read_csv(r'L:\PD\Experiment\Maturation_37\Output_Data\CSV\raw_data_normed_median_dmso_snca.csv')

# # Rename some data tags for consistency (if necessary)
# df.loc[df['tags'].str.contains('Prostratin;WT'), 'tags'] = 'WT+Pro'  

# Overview of data categories
df_columns = df.columns

# Plot info
data=df
x_grouping ="tags"
datacolumn = "sum_intensity_TH_SNCA_SNCA"
plot_order=("Ctrl;DMSO", "Ctrl;PEP005", "Ctrl;Prostratin", "Tripli;DMSO", "Tripli;PEP005", "Tripli;Prostratin")

# Load the above data into `dabest`.
df_dabest = dabest.load(data=data, x=x_grouping, y=datacolumn,
                          idx=plot_order)

# Produce a Cumming estimation plot. Choose "median" instead of "mean" for skewed distributions. 
plt.style.use("seaborn-ticks")

f = df_dabest.median_diff.plot(raw_marker_size=3); #contrast_ylim=(-0.15, 0.4)

rawswarm_axes = f.axes[0]
contrast_axes = f.axes[1]

rawswarm_axes.set_xticklabels(rawswarm_axes.get_xticklabels(),rotation=45, ha='right', rotation_mode="anchor")
contrast_axes.set_xticklabels(contrast_axes.get_xticklabels(),rotation=45, ha='right', rotation_mode="anchor")


"""
The tutorial can be found here:
    
    https://acclab.github.io/DABEST-python-docs/tutorial.html
    
"""