"""
Johannes Wilbertz, PhD (Ksilink) 30 Sep 2021

- Script to median normalize a dataframe based on a reference group (i.e. WT DMSO)
- Each plate within the dataframe is treated separetly

"""

import numpy as np
import pandas as pd
import os

# Merge all relevant .csv files from all subdirectories into a single dataframe
path = "B:\Project - PD\LRRK2 paper\RawData_7.csv"
output_path = r'L:\PROJECTS\PD\Experiment\Maturation_61\NormalizedData_Exp61-63\220607_PD_Exp61-63_normalized.csv'
appended_data = []

df  = pd.read_csv(path)

print(df['tags'].unique())
print(df['Plate'].unique())

# To which reference group should everything be normalized?
reference_group = "LRRK2 WT: DoD 37"

# Get titles of numeric columns for normalization
df_numbers = df.select_dtypes(include=[np.number])
df_numbers_col_names = df_numbers.columns.values.tolist()

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