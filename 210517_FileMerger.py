"""
Merge and sace to csv files according to LDA notbook requirements.

17 May 2021

J. Wilbertz, Ksilink
"""
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Generate dataframes containing individual replicates
df1 = pd.read_csv(r'L:\PD\Experiment\Maturation_35\Output_Data_NewFeatures\210706_NewFeatures_2\agSNCA-Mat35-x40-TH-SNCA-map2-p1-AV.csv')
df2 = pd.read_csv(r'L:\PD\Experiment\Maturation_35\Output_Data_NewFeatures\210706_NewFeatures_2\agSNCA-Mat35-x40-TH-SNCA-map2-p2-AV.csv')
df3 = pd.read_csv(r'L:\PD\Experiment\Maturation_35\Output_Data_NewFeatures\210706_NewFeatures_2\agSNCA-Mat35-x40-TH-SNCA-map2-p3-AV.csv')

# # Add plate/replicate number
# df1.loc[df1['Plate'].str.contains('p2'), 'Plate'] = 'p2'
# df2.loc[df2['Plate'].str.contains('p3'), 'Plate'] = 'p3'

# Append dataframes and adjust organisation for plotting
df = df1.append([df2, df3])

df["tags"]= df["tags"].replace('Prostratin;WT', 'WT + Pro')
df["tags"]= df["tags"].replace('PEP005;WT', 'WT + PEP')
df["tags"]= df["tags"].replace('Prostratin;Mut', 'G2019S + Pro')
df["tags"]= df["tags"].replace('Mut;Prostratin', 'G2019S + Pro')
df["tags"]= df["tags"].replace('Mut;PEP005', 'G2019S + PEP')
df["tags"]= df["tags"].replace('PEP005;Mut', 'G2019S + PEP')
df["tags"]= df["tags"].replace('Mut', 'G2019S + DMSO')
df["tags"]= df["tags"].replace('WT', 'WT + DMSO')

df.drop(['SNCA#20210706_184712', 'SNCA#20210706_184714', 'SNCA#20210706_184713'], axis=1, inplace=True)

# Save as CSV file
df.to_csv(r'L:\PD\Experiment\Maturation_35\Output_Data_NewFeatures\210706_NewFeatures_2\Combined_p1_p2_p3.csv', index = False)