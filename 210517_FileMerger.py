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
df1 = pd.read_csv(r'C:\Users\JohannesWilbertz\Desktop\CheckoutResults\SNCA\Checkout_Results\Default\agSNCA-Mat36-x40-TH-SNCA-map2-p2-AV.csv')
df2 = pd.read_csv(r'C:\Users\JohannesWilbertz\Desktop\CheckoutResults\SNCA\Checkout_Results\Default\agSNCA-Mat36-x40-TH-SNCA-map2-p3-AV.csv')

# Add plate/replicate number
df1.loc[df1['Plate'].str.contains('p2'), 'Plate'] = 'p2'
df2.loc[df2['Plate'].str.contains('p3'), 'Plate'] = 'p3'

# Append dataframes and adjust organisation for plotting
combined = df1.append(df2)

# Add required columns with Metadata from tags
combined.loc[combined['tags'].str.fullmatch('Tripli'), 'tags'] = 'Mut;DMSO'
combined.loc[combined['tags'].str.fullmatch('Ctrl'), 'tags'] = 'WT;DMSO'

combined.loc[combined['tags'].str.contains('Prostratin;Tripli'), 'tags'] = 'Mut;Prostratin'
combined.loc[combined['tags'].str.contains('PEP005;Tripli'), 'tags'] = 'Mut;PEP005'
combined.loc[combined['tags'].str.contains('Ctrl;Prostratin'), 'tags'] = 'WT;Prostratin'
combined.loc[combined['tags'].str.contains('Ctrl;PEP005'), 'tags'] = 'WT;PEP005'

# combined.loc[combined['tags'].str.contains('PEP005'), 'CPD_ID'] = 'PEP005'
# combined.loc[combined['tags'].str.contains('Prostratin'), 'CPD_ID'] = 'Prostratin'
# combined.loc[combined['tags'].str.contains('DMSO'), 'CPD_ID'] = 'DMSO'
               
# combined['Experiment'] = combined['Plate'] + '_' + combined['Class']
# combined['Exp_CPD'] = combined['Experiment'] + '_' + combined['CPD_ID']

# Save as CSV file
combined.to_csv(r'C:\Users\JohannesWilbertz\Desktop\CheckoutResults\SNCA\Checkout_Results\Default\Combined_p2_p3.csv', index = False)