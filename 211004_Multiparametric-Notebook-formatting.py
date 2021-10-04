"""
Merge and save Checkout raw data output to .csv file according to 
multiparametric (LDA) notbook requirements. Normalization will be done in
notebook.

4 October 2021

J. Wilbertz, Ksilink
"""
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Generate dataframes containing individual replicates
df1 = pd.read_csv(r'L:\PD\Experiment\Maturation_44\210929_Exp44_p1\agLRRK2-Mat44-x40-TH-SNCA-MAP2-p1-AV.csv')
df2 = pd.read_csv(r'L:\PD\Experiment\Maturation_44\210929_Exp44_p2\agLRRK2-Mat44-x40-TH-SNCA-MAP2-p2-AV.csv')

df = df1.append(df2)

df.drop(['SNCA#20210929_130725', 'SNCA#20210929_160423'], axis=1, inplace=True)

df['tags'] = df['tags'].astype(str) + ';DMSO'

df_columns = df.columns

df_tags = df['tags'].unique().tolist()

conditions, cpds = zip(*(df['tags'].str.split(';')))

list_class =[]
list_cpds = []

for i, v in enumerate(conditions):
	if cpds[i] in ['Ctrl c89bmS4', 'LRRK2 WT', 'LRRK2 Mut', 'SNCA WT ', 'Ctrl Gibco']:
			list_class.append(cpds[i])
	if cpds[i] in ['DMSO']:
			list_cpds.append(cpds[i])
	if conditions[i] in ['Ctrl c89bmS4', 'LRRK2 WT', 'LRRK2 Mut', 'SNCA WT ', 'Ctrl Gibco']:
			list_class.append(conditions[i])
	if conditions[i] in ['DMSO']:
			list_cpds.append(conditions[i])

conditions = list_class
cpds = list_cpds

df['Class'] = conditions
df['CPD_ID'] = cpds

list_pl = []
list_exp = []
list_exp_cpd = []

plates = df['Plate']

for i, pl in enumerate(plates):
    my_plate = ''
    if 'p1' in pl:
        my_plate = 'p1'
    if 'p2' in pl:
        my_plate = 'p2'

    list_pl.append(my_plate)
    list_exp.append(my_plate+'_'+conditions[i])
    list_exp_cpd.append(my_plate+'_'+conditions[i]+'_'+cpds[i])

df['Plate'] = list_pl
df['Experiment'] = list_exp
df['Exp_CPD'] = list_exp_cpd

# print(df['Plate'].unique())
# print(df['Experiment'].unique())
# print(df['Exp_CPD'].unique())

df_numbers = df.select_dtypes(include=[np.number])
df_col_names = df_numbers.columns.tolist()

print(df_col_names)

# df = df.fillna(0)

# Save as CSV file
df.to_csv(r'L:TEST\PD_Exp44_p1_p2_LDA-ready.csv', index = False)
