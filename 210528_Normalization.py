"""
Normalization per column according to data category of choice (i.e. WT, DMSO)

Type of normalization can be chosen.

Created by Loic Cousin, PhD @ Ksilink and modified by Johannes Wilbertz, PhD @ Ksilink

"""

import pandas as pd
import numpy as np

def norm_centered_by_plate(df, norm_type, neg, pos):
	df['Plate2'] = df.index.get_level_values('Plate')

	print(df['tags'].unique())
	print(df['Plate2'].unique())

	list_df_group = []

	for plate in np.unique(df['Plate2']):
		df_plate = df[df['Plate2'] == plate]

		WT_group = pd.DataFrame(df_plate[df_plate['tags'] == neg])
		GS_group = pd.DataFrame(df_plate[df_plate['tags'] == pos])

		for column in df_plate:
			if column != 'Plate2' and column != 'tags':

				mean_wt = np.mean(WT_group[column])
				mean_gs = np.mean(GS_group[column])
				median_wt = np.median(WT_group[column])

				std_wt = np.std(WT_group[column])
				std_gs = np.std(GS_group[column])

				min_5sigma_wt = 0.
				min_5sigma_gs = 0.
				max_5sigma_wt = 0.
				max_5sigma_gs = 0.

				if norm_type == 'std':
					if std_wt == 0:
						std_wt = 1.

					print(column, mean_wt, std_wt)

					df_plate[column] = (df_plate[column] - median_wt) / std_wt

				if norm_type == 'std_quartiles':

					diff_quartile_wt = np.abs(np.quantile(WT_group[column], 0.75) - np.quantile(WT_group[column], 0.25))

					if diff_quartile_wt == 0:
						diff_quartile_wt = 1.

					print(column, median_wt, diff_quartile_wt)

					df_plate[column] = (df_plate[column] - median_wt) / diff_quartile_wt

				if norm_type == 'mean':
					if mean_wt == 0:
						mean_wt = 1

					print(column, mean_wt)
					df_plate[column] = df_plate[column] / mean_wt

				if norm_type == 'median':
					if median_wt == 0:
						median_wt = 1

					print(column, median_wt)
					df_plate[column] = df_plate[column] / median_wt

				if norm_type == 'min_max':

					min_wt = np.min(WT_group[column])
					max_wt = np.max(WT_group[column])

					min_gs = np.min(GS_group[column])
					max_gs = np.max(GS_group[column])

					sigma = 5

					if mean_wt < mean_gs:
						min_5sigma_wt = mean_wt - sigma * std_wt
						max_5sigma_gs = mean_gs + sigma * std_gs

						if min_5sigma_wt > min_wt:
							min_wt = min_5sigma_wt

						if max_5sigma_gs < max_gs:
							max_gs = max_5sigma_gs

					if mean_gs < mean_wt:
						min_5sigma_gs = mean_gs - sigma * std_gs
						max_5sigma_wt = mean_wt + sigma * std_wt

						if min_5sigma_gs > min_gs:
							min_gs = min_5sigma_gs

						if max_5sigma_wt < max_wt:
							max_wt = max_5sigma_wt

					df_plate[column] -= min_wt

					if (max_wt - min_wt) != 0:
						df_plate[column] /= (max_wt - min_wt)
					else:
						df_plate[column] /= 1.

		list_df_group.append(df_plate)

	normalized_df = pd.concat(list_df_group)
	normalized_df.drop('Plate2', axis=1, inplace=True)

	print('[DATA NORMALIZED]')

	return normalized_df


def remove_outliers_by_plate(df_in, sigma, method, neg, pos):
	print('[DATA BEFORE OUTLIERS:]', len(df_in))

	# numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
	# des = df_in.select_dtypes(include=numerics).columns

	list_df_plate = []

	df_in['Plate2'] = df_in.index.get_level_values('Plate')
	for plate in np.unique(df_in['Plate2']):

		df_plate = pd.DataFrame(df_in[df_in['Plate2'] == plate])

		w_t = df_plate[df_plate['tags'] == neg]
		g_s = df_plate[df_plate['tags'] == pos]
		cpd = df_plate[~((df_plate['tags'] == neg) | (df_plate['tags'] == pos))]

		print('CPDS LEN = ', len(cpd))

		w_t = w_t.drop(['Plate2', 'tags'], axis=1, errors='ignore')
		g_s = g_s.drop(['Plate2', 'tags'], axis=1, errors='ignore')
		cpd = cpd.drop(['Plate2'], axis=1, errors='ignore')

		w_t = w_t.astype(np.float64)
		g_s = g_s.astype(np.float64)

		# print(w_t)

		if method == 'std':
			w_t = w_t[np.abs(w_t - w_t.mean()) / w_t.std() <= sigma].dropna()
			g_s = g_s[np.abs(g_s - g_s.mean()) / g_s.std() <= sigma].dropna()
		elif method == 'std_quartiles':
			w_t = w_t[np.abs(w_t - w_t.median()) <= sigma * np.abs(w_t.quantile(0.75) - w_t.quantile(0.25))].dropna()
			g_s = g_s[np.abs(g_s - g_s.median()) <= sigma * np.abs(g_s.quantile(0.75) - g_s.quantile(0.25))].dropna()

		w_t['tags'] = neg
		g_s['tags'] = pos

		df = pd.concat([w_t, g_s, cpd])
		list_df_plate.append(df)

	df_concat = pd.concat(list_df_plate)

	print('[DATA AFTER OUTLIERS:]', df_concat.shape[0])

	return df_concat

if __name__ == '__main__':

	df = pd.read_csv(r"L:\PD\Experiment\Maturation_35\Output_Data_NewFeatures\210706_NewFeatures_2\Combined_p1_p2_p3.csv", index_col=['Plate', 'Well'])

   # df.drop(['SNCA#20210706_184712', 'SNCA#20210706_184714', 'SNCA#20210706_184713'], axis=1, inplace=True)

	sigma = 3
	# method = "std"
	# method = "mean"
	method = "median"

	# df_wo_outliers = remove_outliers_by_plate(df, sigma=sigma, method=method, neg="A18944;no treat", pos="OtherClass")

	print('Sep')

	df_norm = norm_centered_by_plate(df, norm_type=method, neg="WT + DMSO", pos="G2019S + DMSO")

	print(df_norm)

# Save as CSV file
df_norm.to_csv(r'L:\PD\Experiment\Maturation_35\Output_Data_NewFeatures\210706_NewFeatures_2\Combined_p1_p2_p3_normed.csv', index = False)