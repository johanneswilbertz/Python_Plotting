"""
Johannes Wilbertz, PhD (Ksilink) 30 Sep 2021

- Script to median normalize a dataframe based on a reference group (i.e. WT DMSO)
- Each plate within the dataframe is treated separetly

"""

import numpy as np
import pandas as pd
import os

# Merge all relevant .csv files from all subdirectories into a single dataframe
path = "L:\\PROJECTS\\PD\\Experiment\\Maturation_61_2\\220530_PD_Exp61\\"
output_path = r'L:\PROJECTS\PD\Experiment\Maturation_61_2\220530_PD_Exp61\220531_PD_Exp61_normalized.csv'
appended_data = []

for filename in os.listdir(path):
    if '.fth' in filename and filename[0:2]=='ag':
        df_current  = pd.read_feather(path+'\\' + '\\'+ filename)
        #print(len(df_current.index))
        #df_current['Directory'] = directory
        appended_data.append(df_current)
df = pd.concat(appended_data).reset_index()
df = df.drop(columns='index')

# To which reference group should everything be normalized?
reference_group = "DMSO CTRL;Tripli"

# Get titles of numeric columns for normalization
df_numbers = df.select_dtypes(include=[np.number])
df_numbers_col_names = df_numbers.columns.values.tolist()

print(df['tags'].unique())
print(df['Plate'].unique())

list_df_group = []

# First loop through all "Plates"
for plate in np.unique(df['Plate']):
	df_plate = df[df['Plate'] == plate]

	ref_group = pd.DataFrame(df_plate[df_plate['tags'] == reference_group])
   
# Second loop through all "Columns" within current "Plate"
	for column in df_plate:
		if column in df_numbers_col_names: # Normalize numeric cols only

			median_ref_group = np.median(ref_group[column]) # Calculate median per column (for reference group only)

			if median_ref_group == 0:
				median_ref_group = 1

			df_plate[column] = df_plate[column] / median_ref_group # Normalize each column with median of reference group

	list_df_group.append(df_plate) 

normalized_df = pd.concat(list_df_group)

# Check if normalization was done correctly based on plate level and reference group. Median of reference group must == 1.
df_grouped = normalized_df.groupby(['Plate', 'tags']).median()

# Save as CSV file
normalized_df.to_csv(output_path, index = False)