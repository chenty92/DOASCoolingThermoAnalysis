"""
@author: Tianyi Chen
Version 1.0, 05/02/2018, Framework completed.
"""

import numpy as np

# ===================================================================================
# Analyze EnergyPlus Results (ouput file .csv) for a Single Zone
# ===================================================================================

class EPlusOutput():

	def __init__(self, csvFilePath):
		"""
		Initialization of a new class object
		@param: csvFilePath: str, add .csv file path to the python environment
		"""
		self.path = csvFilePath

	def readFile(self):
		"""
		Read .csv file and prepare the data for post-processing
		@return:
		time: list[str], len(time) = 8760, format as ' 01/01  01:00:00'
		T_out: list[float], len(lighting_raw) = 8760, hourly lights electric energy (J)
		"""
		with open(self.path) as f:
			# Prepare the empty list
			time = []
			T_out = []
			RH_out = []
			T_in = []
			RH_in = []
			sensibleCooling = []
			latentCooling = []
			cooling = []
			# Find the index
			line = f.readline()
			title = line.split(',')
			idx_time = title.index('Date/Time')
			idx_T_out = title.index('Environment:Site Outdoor Air Drybulb Temperature [C](Hourly)')
			idx_RH_out = title.index('Environment:Site Outdoor Air Relative Humidity [%](Hourly)')
			idx_T_in = title.index('LOE:Zone Air Temperature [C](Hourly)')
			idx_RH_in = title.index('LOE:Zone Air Relative Humidity [%](Hourly)')			
			idx_cooling = title.index('LOE IDEAL LOADS AIR:Zone Ideal Loads Supply Air Total Cooling Energy [J](Hourly)')
			idx_sensibleCooling = title.index('LOE IDEAL LOADS AIR:Zone Ideal Loads Supply Air Sensible Cooling Energy [J](Hourly)')
			idx_latentCooling = title.index('LOE IDEAL LOADS AIR:Zone Ideal Loads Supply Air Latent Cooling Energy [J](Hourly)')
			# Collect hourly data         
			for line in f:
				catalog = line.split(',')
				time.append(catalog[idx_time])
				T_out.append(float(catalog[idx_T_out]))
				RH_out.append(float(catalog[idx_RH_out]))
				T_in.append(float(catalog[idx_T_in]))
				RH_in.append(float(catalog[idx_RH_in]))
				cooling.append(float(catalog[idx_cooling]))
				sensibleCooling.append(float(catalog[idx_sensibleCooling]))
				latentCooling.append(float(catalog[idx_latentCooling]))
				
		return time, T_out, RH_out, T_in, RH_in, cooling, sensibleCooling, latentCooling