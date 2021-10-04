"""
Johannes Wilbertz, PhD (Ksilink) 30 Sep 2021

- Script to median normalize a dataframe based on a reference group (i.e. WT DMSO)
- Each plate within the dataframe is treated separetly

"""

import numpy as np
import pandas as pd

# Generate dataframe
df = pd.read_csv(r'L:\PD\Experiment\Maturation_46\PD_Exp46_p1_p2_MT-TMRM_LDA-ready.csv')

# To which reference group should everything be normalized?
reference_group = "Ctrl;DMSO"

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
normalized_df.to_csv(r'L:\PD\Experiment\Maturation_46\PD_Exp46_p1_p2_normalized.csv', index = False)