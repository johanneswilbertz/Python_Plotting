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

# Generate dataframe
df = pd.read_csv(r'L:\PD\Experiment\Maturation_26\Output_Data\CSV\raw_data_normed.csv') #_normed

# Data wrangling and renaming
df = df[df.tags.str.contains("no treat")] # Add ~df.tags.str.contains("no treat") to exclude, remove to include
df = df[~df.tags.str.contains("EDI001A1/2")]
df.loc[df['tags'].str.contains('A18944'), 'tags'] = 'Control 1'
df.loc[df['tags'].str.contains('c89bmS4'), 'tags'] = 'Control 2'
df.loc[df['tags'].str.contains('Edi001A'), 'tags'] = 'Tripli'

# Overview of data categories
df_columns = df.columns

# Load the above data into `dabest`.
df_dabest = dabest.load(data=df, x="tags", y="surface_p-SNCA",
                          idx=("Control 1", "Control 2", "Tripli"))

# Produce a Cumming estimation plot.
df_dabest.mean_diff.plot();

"""
The tutorial can be found here:
    
    https://acclab.github.io/DABEST-python-docs/tutorial.html
    
"""