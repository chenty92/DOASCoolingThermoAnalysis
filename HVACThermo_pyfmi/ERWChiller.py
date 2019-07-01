"""
@author: Tianyi Chen
Version 1.0, 12/04/2018, Framework completed.
"""

from pyfmi import load_fmu
import numpy as np
import EPlusOutput as epo

# =============================================================================
# Load Modelica Model and Run Parametric Study Cases
# =============================================================================

class runERWChillerModel():

	def __init__(self, csvFilePath, modelPath, idx_start, idx_end, V_a, V_w, epsilon, epsilon_s, epsilon_l, r, 
		phi_min, t_lcmin, deltaCOP, deltaT, deltaw, deltaP, deltaQ, deltaWc, W_fan, P_fan):
		"""
		Inputs:
		# 	csvFilePath: str, add .csv file path to the python environment
		# 	modelPath: str, add .fmu file path to the python environment
		# 	idx_start: int, index of start time for simulation
		# 	idx_end: int, index of end time for simulation
		#   V_a: float, volumetric flow rate of moist air, unit: m3/s
		#   V_w: float, volumetric flow rate of water, unit: m3/s
		#   epsilon: float, performance of HX in chiller
		#   epsilon_s: float, performance of ERW in sensible heat recovery
		#   epsilon_l: float, performance of ERW in latent heat recovery
		#	r: float, ratio of mass flow rate of exhaust air in ERW to supply air
		#	unit: float, number of number of parallel terminals
		#	phi_min: float, minimum absolute humidity for supply air
		#	t_lcmin: float, minimum leaving chilled water temperature, unit: C
		#	deltaCOP: float, threshold control on chiller
		#	deltaT: float, temperature threshold control on/off for the ERV system, unit: K
		#	deltaw: float, humidity threshold control on/off for the ERV system
		#	deltaP: float, pressure threshold control on/off for the ERV system, unit: Pa
		#	deltaQ: float, threshold control on total cooling load, unit: W
		#	deltaWc: float, threshold control on minimum cooling work, unit: W
		#	W_fan: float, fan power, unit: W	
		"""
		self.path = csvFilePath
		self.model = modelPath
		self.start = idx_start
		self.end = idx_end
		self.V_a = V_a
		self.V_w = V_w
		self.epsilon = epsilon
		self.epsilon_s = epsilon_s
		self.epsilon_l = epsilon_l
		self.phi_min = phi_min
		self.t_lcmin = t_lcmin
		self.r = r
		self.deltaT = deltaT
		self.deltaw = deltaw
		self.deltaP = deltaP
		self.deltaCOP = deltaCOP
		self.deltaQ = deltaQ
		self.deltaWc = deltaWc
		self.W_fan = W_fan
		self.P_fan = P_fan

	def runModelicaModel(self): 
		"""
		Outputs:
		#   COP: list[float], time series solutions of COP for the baseline chiller system in kW, len(COP) = idx_end - idx_start
		#   T_lc: list[float], time series solutions of leaving chilled water temperature in C, len(T_lc) = idx_end - idx_start
		#   T_ec: list[float], time series solutions of entering chilled water temperature in C, len(T_ec) = idx_end - idx_start
		#   Q_c: list[float], time series solutions of cooling load for the baseline chiller system in kW, len(Q_c) = idx_end - idx_start
		"""

		# Read the EPlus output file
		data = epo.EPlusOutput(self.path)
		time, T_out, RH_out, T_in, RH_in, cooling, sensibleCooling, latentCooling = data.readFile()

		# Prepare the result output
		COP = [0] * (self.end - self.start)
		T_lc = [0] * (self.end - self.start)
		T_ec = [0] * (self.end - self.start)
		Q_c = [0] * (self.end - self.start)

		# Load the model
		model = load_fmu(self.model)

		# Run quasi-steady model for each time

		for i in range(self.start, self.end):

			# Prepare the inputs
			var = ['t_o', 'rh_o', 't_i', 'rh_i', 'q_lat']
			t = [0]
			t_o = [T_out[i]]
			rh_o = [RH_out[i]]
			t_i = [T_in[i]]
			rh_i = [RH_in[i]]
			q_lat = [latentCooling[i]]

			# Reset model inputs and parameters
			model.reset()

			# Model parameter setting
			epsilon_s = self.epsilon_s
			model.set('epsilon_s', epsilon_s)
			epsilon_l = self.epsilon_l
			model.set('epsilon_l', epsilon_l)
			V_a = self.V_a
			model.set('V_a', V_a)
			V_w = self.V_w
			model.set('V_w', V_w)
			epsilon = self.epsilon
			model.set('epsilon', epsilon)
			phi_min = self.phi_min
			model.set('phi_min', phi_min)
			t_lcmin = self.t_lcmin
			model.set('t_lcmin', t_lcmin)
			deltaCOP = self.deltaCOP
			model.set('deltaCOP', deltaCOP)
			r = self.r
			model.set('r', r)
			deltaT = self.deltaT
			model.set('deltaT', deltaT)
			deltaw = self.deltaw
			model.set('deltaT', deltaw)
			deltaP = self.deltaP
			model.set('deltaP', deltaP)
			deltaQ = self.deltaQ
			model.set('deltaQ', deltaQ)
			deltaWc = self.deltaWc
			model.set('deltaWc', deltaWc)
			W_fan = self.W_fan
			model.set('W_fan', W_fan)
			P_fan = self.P_fan
			model.set('P_fan', P_fan)

			# Update the inputs
			inputArr = np.transpose(np.array([t, t_o, rh_o, t_i, rh_i, q_lat]))
			input_object = (var, inputArr)

			# Update the model setting
			simulate_opts = model.simulate_options()
			simulate_opts['CVode_options']['rtol'] = 1e-6

			# Simulate the model
			try:
				res = model.simulate(start_time = 0, final_time = 1, input = input_object, options = simulate_opts)

				# Results
				COP[i - self.start] = res['COP'][-1]
				Q_c[i - self.start] = res['q_c'][-1]

			except:
				print('----------Warning: Iteration ' + str(i) + ' has initialization error!----------')
				COP[i - self.start] = None
				Q_c[i - self.start] = None	

		return COP, Q_c