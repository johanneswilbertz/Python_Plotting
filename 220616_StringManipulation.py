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
from sys import exit

# Load data
path = r'L:\PROJECTS\PD\Experiment\Maturation_61_2\220530_PD_Exp61\220531_PD_Exp61.csv'
df  = pd.read_csv(path)

# Get unique data tags/categories/columns
df_columns = df.columns
unique_tags = df['tags'].unique()
print("ORIGINAL: ", unique_tags)

# Correct tags
df['tags'] = df['tags'].str.replace('rawA', 'Row A')
df['tags'] = df['tags'].str.replace('rawP', 'Row P')
df['tags'] = df['tags'].str.replace('raw A', 'Row A')
df['tags'] = df['tags'].str.replace('raw P', 'Row P')

# Remove all row annotations
df['tags'] = df['tags'].str.replace(" Row A", "")
df['tags'] = df['tags'].str.replace(" Row P", "")

# Simplify plate names
df['Plate'] = df['Plate'].str.replace('Mat61 Tripli DMSO','Plate 1')
df['Plate'] = df['Plate'].str.replace('Mat61 Tripli Prostratin','Plate 2')
df['Plate'] = df['Plate'].str.replace('Mat61-WT','Plate 3')

# Extract rows and columns
df['Row'] = df['Well'].astype(str).str[0]
df['Column'] = df['Well'].astype(str).str[1:3]

# Delete all WT wells
df = df[df['tags'].str.contains('WT')==False]

# Add row/column info to tag based on specific rows or columns
df['tags'] = np.where(df['Row'] == 'A', df['tags']+" "+ 'Row A', df['tags'])
df['tags'] = np.where(df['Row'] == 'P', df['tags']+" "+ 'Row P', df['tags'])
df['tags'] = np.where(df['Row'] == 'B', df['tags']+" "+ 'Row B', df['tags'])
df['tags'] = np.where(df['Row'] == 'O', df['tags']+" "+ 'Row O', df['tags'])
df['tags'] = np.where(df['Column'] == '02', df['tags']+" "+ 'Column 02', df['tags'])
df['tags'] = np.where(df['Column'] == '23', df['tags']+" "+ 'Column 23', df['tags'])

unique_tags = df['tags'].unique()
print("EDITED: ", unique_tags)
print(df['tags'].value_counts())

#output_path = r'L:\PROJECTS\PD\Experiment\Maturation_61_2\220530_PD_Exp61\220531_PD_Exp61_IntraPlateModified.csv'
#df.to_csv(output_path, index = False)