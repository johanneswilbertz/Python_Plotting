# -*- coding: utf-8 -*-
"""
Johannes Wilbertz @ Ksilink, 01 June 2021

Based on:

Moving beyond P values: Everyday data analysis with estimation plots
Joses Ho, Tayfun Tumkaya, Sameer Aryal, Hyungwon Choi, Adam Claridge-Chang
Nature Methods 2019, 1548-7105. 10.1038/s41592-019-0470-3

"""
import pandas as pd
import numpy as np
from scipy import stats
import dabest
import matplotlib.pyplot as plt

# Generate dataframes containing individual replicates
df = pd.read_csv(r'L:\PD\Experiment\Maturation_25\Output_Data\CSV\raw_data_normed_median_dmso.csv')

# Rename some data tags for consistency (if necessary)
df["tags"]= df["tags"].replace('Mut;PFE', 'PFE;Mut') 

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

# Which columns are left?
df_columns = df.columns

# Plot info
data=df
x_grouping ="tags"
datacolumn = "sum_intensity_per_nuclei_SNCA"
plot_order=("DMSO;WT", "GNE;WT", "PFE;WT", "MLi;WT", "DMSO;Mut", "GNE;Mut", "PFE;Mut", "MLi;Mut")

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