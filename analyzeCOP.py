"""
@author: Tianyi Chen
Version 1.0, 12/27/2018, Framework completed.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

##### General Settings
# Ventilation rate
TFA = 1750						# Unit: m2
V_a = 0.55 * 1e-3 * TFA 		# Unit: m3/s
# Supply/exhaust fan
W_fan = 1.5 * V_a * 1e3			# Unit: W

##### Load file
filePath = 'D:/tianyich/Dropbox (MIT)/PhD Work/PassiveStrategies4Cooling/'
fileName = 'Mumbai_Hybrid_HighTM'

file = filePath + fileName + '_res.csv'

# file = 'boston_loE_res.csv'
df = pd.read_csv(file, header = 0)

##### Calculate Monthly working hours
numHrs = []
numWkgHrs = []
Q_c = []
time = df['Date/Time'].tolist()
# List all the non-zero entry row index
idx = list(df['Q_c (kW)'] != 0)
# Iterative thru the rows
for i in range(1, 13):
	idx_Q_c = [False] * df.shape[0]
	if i == 1:
		idx_start = 0
	else:
		idx_start = time.index(' %02d/01  01:00:00' % i)
	if i == 12:
		idx_end = -1
		# numHrsInMonth = float(df.index.max() - idx_start)
		idx_Q_c[idx_start: df.index.max()] = idx[idx_start: df.index.max()]		
	else:
		idx_end = time.index(' %02d/01  01:00:00' % (i + 1))
		# numHrsInMonth = float(idx_end - idx_start)
		idx_Q_c[idx_start: idx_end] = idx[idx_start: idx_end]
	
	numHrsInMonth = float((df['Q_lat (kW)'][idx_start: idx_end] != 0).sum())
	numHrs.append(numHrsInMonth)
	numWkgHrsinMonth = float(((df['Q_c (kW)'][idx_start: idx_end] != 0) & 
							  (df['Q_lat (kW)'][idx_start: idx_end] != 0)).sum())
	# numWkgHrs.append(numHrsInMonth - (df['Q_c (kW)'][idx_start: idx_end] == 0).sum())	
	numWkgHrs.append(numWkgHrsinMonth)
	idx_Q_c[idx_start: idx_end] = ((df['Q_c (kW)'][idx_start: idx_end] != 0) & 
							  	   (df['Q_lat (kW)'][idx_start: idx_end] != 0))
	Q_c.append(df['Q_c (kW)'][idx_Q_c].sum()/TFA)
fracWkgHrs = np.array(numWkgHrs) / np.array(numHrs)
print fracWkgHrs, numHrs, Q_c


##### Remove all the rows with Q_c = 0
df = df[df['Q_c (kW)'] != 0]
df = df[df['Q_lat (kW)'] != 0]

##### Calculate the work input for each HVAC Cooling system
df['w_bcs (kW)'] = df['Q_c (kW)'] / df['COP_bcs']
df['w_dcs (kW)'] = df['Q_c (kW)'] / df['COP_dcs']
df['w_mcs (kW)'] = df['Q_c (kW)'] / df['COP_mcs']
df['w_mcms (kW)'] = df['Q_c (kW)'] / df['COP_mcms']
df['w_ecs (kW)'] = df['Q_c (kW)'] / df['COP_ecs']
df['w_hcs (kW)'] = df['Q_c (kW)'] / df['COP_hcs']


##### Calculate the annual average COP for each system
Q_c_tot = sum(df['Q_c (kW)'])   # Unit: kWh
COP_avg_bcs = sum(df['Q_c (kW)']) / sum(df['w_bcs (kW)'])
COP_avg_dcs = sum(df['Q_c (kW)']) / sum(df['w_dcs (kW)'])
COP_avg_mcs = sum(df['Q_c (kW)']) / sum(df['w_mcs (kW)'])
COP_avg_mcms = sum(df['Q_c (kW)']) / sum(df['w_mcms (kW)'])
COP_avg_ecs = sum(df['Q_c (kW)']) / sum(df['w_ecs (kW)'])
COP_avg_hcs = sum(df['Q_c (kW)']) / sum(df['w_hcs (kW)'])

print df.shape[0]
print (Q_c_tot, COP_avg_bcs, COP_avg_dcs, COP_avg_mcs, COP_avg_mcms, COP_avg_ecs, COP_avg_hcs)
# Bos: (15000.96, 2.73, 5.50, 5.25, 5.37, 2.84, 2.70)
# Bos_Summer: (12274.13, 2.83, 5.73, 5.69, 5.71, 2.96, 2.79)
# Mumbai: (185395.34, 3.99, 8.93, 9.11, 10.17, 5.28, 4.03)

# Bos_loE: (5401.88, 3.73, 7.23, 7.99, 6.58, 3.93, 3.74)


##### Remove the fan power of the system and look at COP for each device
def COP_dev(devName, df, W_fan, minWork = 1.0, maxCOP = 50.0):
	"""
	Inputs:
	devName: string, e.g. 'w_bcs (kW)'
	df: dataframe
	minWork: float, default = 0.1 kW
	W_fan: float, unit: kW

	Output:
	COP: dataframe
	"""
	# First filter thru the data, remove cases when sensible cooling dominant
	# (net device work <= minWork)
	d = df[df[devName] >= W_fan + minWork]
	COP = (d['Q_c (kW)'] / (d[devName] - W_fan)).clip(0.0, maxCOP)

	return COP

COP_bcs = COP_dev('w_bcs (kW)', df, W_fan / 1e3)
COP_dcs = COP_dev('w_dcs (kW)', df, W_fan / 1e3)
COP_mcs = COP_dev('w_mcs (kW)', df, W_fan / 1e3)
COP_mcms = COP_dev('w_mcms (kW)', df, W_fan / 1e3)
COP_ecs = COP_dev('w_ecs (kW)', df, W_fan / 1e3)
COP_hcs = COP_dev('w_hcs (kW)', df, W_fan / 1e3)


#########################################################################
##### Plot
# # Plot fraction of working hrs per month
# fig = plt.figure()
# plt.rcParams["font.family"] = "Times New Roman"
# x = range(len(numHrs))
# monthName = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']
# width = 0.5
# plt.bar(x, fracWkgHrs, width, color='b')
# plt.title('DOAS Cooling System Monthly Working Hours')
# plt.xticks([num for num in x], monthName)
# plt.xlabel('Time (Month)')
# plt.ylabel('Fraction of Working Hours per Month')
# fig.savefig('DOAS Cooling System Monthly Working Hours_Boston.png', bbox_inches='tight', dpi=1200)

# # Plot DOAS Cooling Load per Month
# fig = plt.figure()
# plt.rcParams["font.family"] = "Times New Roman"
# x = range(len(numHrs))
# monthName = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']
# width = 0.5
# plt.bar(x, Q_c, width, color='b')
# plt.title('DOAS Cooling System Monthly Working Load')
# plt.xticks([num for num in x], monthName)
# plt.xlabel('Time (Month)')
# plt.ylabel('Working Load per Month ($kWh/m^2$)')
# fig.savefig('DOAS Cooling System Monthly Working Load_Boston.png', bbox_inches='tight', dpi=1200)

# # Plot density of COP distribution
# fig = plt.figure()
# plt.rcParams["font.family"] = "Times New Roman"
# maxRange = max(df['COP_bcs'].max(), df['COP_dcs'].max(), df['COP_mcs'].max(), df['COP_mcms'].max(), df['COP_ecs'].max(), df['COP_hcs'].max())
# sns.distplot(list(df['COP_bcs']), hist=False, kde=True, kde_kws={'color': 'k', 'lw': 1, 'label': 'BCS'}).set(xlim=(0, maxRange)) 
# sns.distplot(list(df['COP_dcs']), hist=False, kde=True, kde_kws={'color': 'r', 'lw': 1, 'label': 'DCS'}).set(xlim=(0, maxRange))
# sns.distplot(list(df['COP_mcs']), hist=False, kde=True, kde_kws={'color': 'b', 'lw': 1, 'label': 'MCCS'}).set(xlim=(0, maxRange))
# sns.distplot(list(df['COP_mcms']), hist=False, kde=True, kde_kws={'color': 'g', 'lw': 1, 'label': 'MCMS'}).set(xlim=(0, maxRange))
# sns.distplot(list(df['COP_ecs']), hist=False, kde=True, kde_kws={'color': 'c', 'lw': 1, 'label': 'ECS'}).set(xlim=(0, maxRange)) 
# sns.distplot(list(df['COP_hcs']), hist=False, kde=True, kde_kws={'color': 'y', 'lw': 1, 'label': 'HCS'}).set(xlim=(0, maxRange))
# plt.title('COP Distribution of HVAC Cooling Systems')
# plt.xlabel('COP')
# plt.ylabel('Density')
# fig.savefig('COP distribution of HVAC Cooling Systems_Boston.png', bbox_inches='tight', dpi=1200)

# # # Plot density of COP distribution
# # fig = plt.figure()
# # plt.rcParams["font.family"] = "Times New Roman"
# # maxRange = max(COP_bcs.max(), COP_dcs.max(), COP_mcs.max(), COP_mcms.max(), COP_ecs.max(), COP_hcs.max())
# # sns.distplot(COP_bcs, hist=False, kde=True, kde_kws={'color': 'k', 'lw': 1, 'label': 'BCS'}).set(xlim=(0, maxRange)) 
# # sns.distplot(COP_dcs, hist=False, kde=True, kde_kws={'color': 'r', 'lw': 1, 'label': 'DCS'}).set(xlim=(0, maxRange))
# # sns.distplot(COP_mcs, hist=False, kde=True, kde_kws={'color': 'b', 'lw': 1, 'label': 'MCCS'}).set(xlim=(0, maxRange))
# # sns.distplot(COP_mcms, hist=False, kde=True, kde_kws={'color': 'g', 'lw': 1, 'label': 'MCMS'}).set(xlim=(0, maxRange))
# # sns.distplot(COP_ecs, hist=False, kde=True, kde_kws={'color': 'c', 'lw': 1, 'label': 'ECS'}).set(xlim=(0, maxRange)) 
# # sns.distplot(COP_hcs, hist=False, kde=True, kde_kws={'color': 'y', 'lw': 1, 'label': 'HCS'}).set(xlim=(0, maxRange))
# # plt.title('COP Distribution of HVAC Cooling Devices')
# # plt.xlabel('COP')
# # plt.ylabel('Density')
# # fig.savefig('COP distribution of HVAC Cooling Devices_Boston.png', bbox_inches='tight', dpi=1200)

###### Plot work distribution
fig = plt.figure()
plt.rcParams["font.family"] = "Times New Roman"
maxRange = 10.0
sns.distplot(df['w_bcs (kW)'], hist=False, kde=True, kde_kws={'color': 'k', 'lw': 1, 'label': 'BCS'}).set(xlim=(0, maxRange)) 
sns.distplot(df['w_dcs (kW)'], hist=False, kde=True, kde_kws={'color': 'r', 'lw': 1, 'label': 'DCS'}).set(xlim=(0, maxRange))
sns.distplot(df['w_mcs (kW)'], hist=False, kde=True, kde_kws={'color': 'b', 'lw': 1, 'label': 'MCCS'}).set(xlim=(0, maxRange))
sns.distplot(df['w_mcms (kW)'], hist=False, kde=True, kde_kws={'color': 'g', 'lw': 1, 'label': 'MCMS'}).set(xlim=(0, maxRange))
sns.distplot(df['w_ecs (kW)'], hist=False, kde=True, kde_kws={'color': 'c', 'lw': 1, 'label': 'ECS'}).set(xlim=(0, maxRange)) 
sns.distplot(df['w_hcs (kW)'], hist=False, kde=True, kde_kws={'color': 'y', 'lw': 1, 'label': 'HCS'}).set(xlim=(0, maxRange))
plt.title('Cooling Work Distribution of HVAC Cooling Systems')
plt.xlabel('Cooling Work (kW)')
plt.ylabel('Density')
fig.savefig('Cooling work distribution of HVAC Cooling Systems_Mumbai.png', bbox_inches='tight', dpi=1200)


##########################################
# # Monthly average COP
# idx_start = 4344
# idx_end = 5087
# idx_end = 6487
# idx_list = [False] * df.shape[0]
# idx_list[idx_start: idx_end] = df['Q_c (kW)'][idx_start: idx_end] != 0
# ##### Remove all the rows with 0
# # df = df[~(df['Q_c (kW)'] == 0).any(axis = 1)]
# f = df[idx_list]

# ##### Calculate the work input for each HVAC Cooling system
# f['w_bcs (kW)'] = f['Q_c (kW)'] / f['COP_bcs']
# f['w_dcs (kW)'] = f['Q_c (kW)'] / f['COP_dcs']
# f['w_mcs (kW)'] = f['Q_c (kW)'] / f['COP_mcs']


# ##### Calculate the annual average COP for each system
# Q_c_tot = sum(f['Q_c (kW)'])   # Unit: kWh
# COP_avg_bcs = sum(f['Q_c (kW)']) / sum(f['w_bcs (kW)'])
# COP_avg_dcs = sum(f['Q_c (kW)']) / sum(f['w_dcs (kW)'])
# COP_avg_mcs = sum(f['Q_c (kW)']) / sum(f['w_mcs (kW)'])

# print (Q_c_tot, COP_avg_bcs, COP_avg_dcs, COP_avg_mcs)
# # Result: (5942.244733132237, 2.204756748634062, 2.6999942516667934, 3.7986251741048838)