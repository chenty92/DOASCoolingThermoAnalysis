"""
@author: Tianyi Chen
Version 1.0, 05/29/2018, Framework completed.
"""

import numpy as np
from pyfmi import load_fmu


# =============================================================================
# Load BaselineChiller Model and Run Parametric Study Cases
# =============================================================================

class runMinWorkModel():

	def __init__(self, T_out, RH_out, T_in, RH_in, latentCooling, modelPath, V_a, phi_min, deltaQ,
		epsilon_s, epsilon_l, deltaP, deltaT, deltaw):
		"""
		Inputs:
		# 	modelPath: str, add .fmu file path to the python environment
		#   V_a: float, volumetric flow rate of moist air, unit: m3/s
		#   V_w: float, volumetric flow rate of water, unit: m3/s
		#   epsilon: float, performance of HX in chiller
		#	unit: float, number of number of parallel terminals
		#	phi_min: float, minimum absolute humidity for supply air
		#	t_lcmin: float, minimum leaving chilled water temperature, unit: C
		#	deltaCOP: float, threshold control on chiller
		#	deltaQ: float, threshold control on total cooling load, unit: W
		#	deltaWc: float, threshold control on minimum cooling work, unit: W
		#	W_fan: float, fan power, unit: W
		"""
		self.T_out = T_out
		self.RH_out = RH_out
		self.T_in = T_in
		self.RH_in = RH_in
		self.q_lat = latentCooling
		self.model = modelPath
		self.V_a = V_a
		self.phi_min = phi_min
		self.deltaQ = deltaQ
		self.epsilon_s = epsilon_s
		self.epsilon_l = epsilon_l
		self.deltaT = deltaT
		self.deltaw = deltaw
		self.deltaP = deltaP

	def runModelicaModel(self): 

		# Prepare the result output
		w_min_sa = np.zeros((len(self.T_out), len(self.RH_out)))
		w_min_l = np.zeros((len(self.T_out), len(self.RH_out)))
		w_min_hw = np.zeros((len(self.T_out), len(self.RH_out)))
		w_min_erw = np.zeros((len(self.T_out), len(self.RH_out)))

		# Load the model
		model = load_fmu(self.model)

		# Run quasi-steady model for each time

		for i in range(len(self.T_out)):
			for j in range(len(self.RH_out)):
				# Prepare the inputs
				var = ['t_o', 'rh_o', 't_i', 'rh_i', 'q_lat']
				t = [0]
				t_o = [self.T_out[i]]
				rh_o = [self.RH_out[j]]
				t_i = [self.T_in[i]]
				rh_i = [self.RH_in[j]]
				q_lat = [self.q_lat[i][j]]

				# Reset model inputs and parameters
				model.reset()

				# Model parameter setting
				V_a = self.V_a
				model.set('V_a', V_a)
				phi_min = self.phi_min
				model.set('phi_min', phi_min)
				deltaQ = self.deltaQ
				model.set('deltaQ', deltaQ)
				epsilon_s = self.epsilon_s
				model.set('epsilon_s', epsilon_s)
				epsilon_l = self.epsilon_l
				model.set('epsilon_l', epsilon_l)
				deltaT = self.deltaT
				model.set('deltaT', deltaT)
				deltaw = self.deltaw
				model.set('deltaw', deltaw)
				deltaP = self.deltaP
				model.set('deltaP', deltaP)

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
					w_min_sa[i][j] = res['w_min_sa'][-1]
					w_min_l[i][j] = res['w_min_l'][-1]
					w_min_hw[i][j] = res['w_min_hw'][-1]
					w_min_erw[i][j] = res['w_min_erw'][-1]
				
				except:
					print('----------Warning: Iteration ' + str(i) + '&' + str(j) + ' has initialization error!----------')
					w_min_sa[i][j] = None
					w_min_l[i][j] = None
					w_min_hw[i][j] = None
					w_min_erw[i][j] = None

		return w_min_sa, w_min_l, w_min_hw, w_min_erw


# =============================================================================
# Load BaselineChiller Model and Run Parametric Study Cases
# =============================================================================

class runBaselineChillerModel():

	def __init__(self, T_out, RH_out, T_in, RH_in, latentCooling, modelPath, V_a, V_w, epsilon, 
		phi_min, t_lcmin, deltaCOP, deltaQ, deltaWc, W_fan):
		"""
		Inputs:
		# 	modelPath: str, add .fmu file path to the python environment
		#   V_a: float, volumetric flow rate of moist air, unit: m3/s
		#   V_w: float, volumetric flow rate of water, unit: m3/s
		#   epsilon: float, performance of HX in chiller
		#	unit: float, number of number of parallel terminals
		#	phi_min: float, minimum absolute humidity for supply air
		#	t_lcmin: float, minimum leaving chilled water temperature, unit: C
		#	deltaCOP: float, threshold control on chiller
		#	deltaQ: float, threshold control on total cooling load, unit: W
		#	deltaWc: float, threshold control on minimum cooling work, unit: W
		#	W_fan: float, fan power, unit: W
		"""
		self.T_out = T_out
		self.RH_out = RH_out
		self.T_in = T_in
		self.RH_in = RH_in
		self.q_lat = latentCooling
		self.model = modelPath
		self.V_a = V_a
		self.V_w = V_w
		self.epsilon = epsilon
		self.phi_min = phi_min
		self.t_lcmin = t_lcmin
		self.deltaCOP = deltaCOP
		self.deltaQ = deltaQ
		self.deltaWc = deltaWc
		self.W_fan = W_fan

	def runModelicaModel(self): 

		# Prepare the result output
		COP = np.zeros((len(self.T_out), len(self.RH_out)))
		Q_c = np.zeros((len(self.T_out), len(self.RH_out)))
		T_lc = np.zeros((len(self.T_out), len(self.RH_out)))
		T_ec = np.zeros((len(self.T_out), len(self.RH_out)))

		# Load the model
		model = load_fmu(self.model)

		# Run quasi-steady model for each time

		for i in range(len(self.T_out)):
			for j in range(len(self.RH_out)):
				# Prepare the inputs
				var = ['t_o', 'rh_o', 't_i', 'rh_i', 'q_lat']
				t = [0]
				t_o = [self.T_out[i]]
				rh_o = [self.RH_out[j]]
				t_i = [self.T_in[i]]
				rh_i = [self.RH_in[j]]
				q_lat = [self.q_lat[i][j]]

				# Reset model inputs and parameters
				model.reset()

				# Model parameter setting
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
				deltaQ = self.deltaQ
				model.set('deltaQ', deltaQ)
				deltaWc = self.deltaWc
				model.set('deltaWc', deltaWc)
				W_fan = self.W_fan
				model.set('W_fan', W_fan)	

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
					COP[i][j] = res['COP'][-1]
					Q_c[i][j] = res['q_c'][-1]

				except:
					print('----------Warning: Iteration ' + str(i) + '&' + str(j) + ' has initialization error!----------')
					COP[i][j] = None
					Q_c[i][j] = None	

		return COP, Q_c



# =============================================================================
# Load HRWChiller Model and Run Parametric Study Cases
# =============================================================================
class runHRWChillerModel():

	def __init__(self, T_out, RH_out, T_in, RH_in, latentCooling, modelPath, V_a, V_w, epsilon, r, 
		phi_min, t_lcmin, deltaCOP, deltaT, deltaQ, deltaWc, W_fan, P_fan):
		"""
		Inputs:
		# 	modelPath: str, add .fmu file path to the python environment
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
		self.T_out = T_out
		self.RH_out = RH_out
		self.T_in = T_in
		self.RH_in = RH_in
		self.q_lat = latentCooling
		self.model = modelPath
		self.V_a = V_a
		self.V_w = V_w
		self.epsilon = epsilon
		self.phi_min = phi_min
		self.t_lcmin = t_lcmin
		self.r = r
		self.deltaT = deltaT
		self.deltaCOP = deltaCOP
		self.deltaQ = deltaQ
		self.deltaWc = deltaWc
		self.W_fan = W_fan
		self.P_fan = P_fan

	def runModelicaModel(self): 

		# Prepare the result output
		COP = np.zeros((len(self.T_out), len(self.RH_out)))
		Q_c = np.zeros((len(self.T_out), len(self.RH_out)))
		T_lc = np.zeros((len(self.T_out), len(self.RH_out)))
		T_ec = np.zeros((len(self.T_out), len(self.RH_out)))

		# Load the model
		model = load_fmu(self.model)

		# Run quasi-steady model for each time

		for i in range(len(self.T_out)):
			for j in range(len(self.RH_out)):
				# Prepare the inputs
				var = ['t_o', 'rh_o', 't_i', 'rh_i', 'q_lat']
				t = [0.0]
				t_o = [self.T_out[i]]
				rh_o = [self.RH_out[j]]
				t_i = [self.T_in[i]]
				rh_i = [self.RH_in[j]]
				q_lat = [self.q_lat[i][j]]

				# Reset model inputs and parameters
				model.reset()

				# Model parameter setting
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
					COP[i][j] = res['COP'][-1]
					Q_c[i][j] = res['q_c'][-1]

				except:
					print('----------Warning: Iteration ' + str(i) + '&' + str(j) + ' has initialization error!----------')
					COP[i][j] = None
					Q_c[i][j] = None	

		return COP, Q_c





# =============================================================================
# Load ERWChiller Model and Run Parametric Study Cases
# =============================================================================
class runERWChillerModel():

	def __init__(self, T_out, RH_out, T_in, RH_in, latentCooling, modelPath, V_a, V_w, epsilon, epsilon_s, epsilon_l, r, 
		phi_min, t_lcmin, deltaCOP, deltaT, deltaw, deltaP, deltaQ, deltaWc, W_fan, P_fan):
		"""
		Inputs:
		# 	modelPath: str, add .fmu file path to the python environment
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
		self.T_out = T_out
		self.RH_out = RH_out
		self.T_in = T_in
		self.RH_in = RH_in
		self.q_lat = latentCooling
		self.model = modelPath
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

		# Prepare the result output
		COP = np.zeros((len(self.T_out), len(self.RH_out)))
		Q_c = np.zeros((len(self.T_out), len(self.RH_out)))
		T_lc = np.zeros((len(self.T_out), len(self.RH_out)))
		T_ec = np.zeros((len(self.T_out), len(self.RH_out)))

		# Load the model
		model = load_fmu(self.model)

		# Run quasi-steady model for each time

		for i in range(len(self.T_out)):
			for j in range(len(self.RH_out)):
				# Prepare the inputs
				var = ['t_o', 'rh_o', 't_i', 'rh_i', 'q_lat']
				t = [0.0]
				t_o = [self.T_out[i]]
				rh_o = [self.RH_out[j]]
				t_i = [self.T_in[i]]
				rh_i = [self.RH_in[j]]
				q_lat = [self.q_lat[i][j]]

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
					COP[i][j] = res['COP'][-1]
					Q_c[i][j] = res['q_c'][-1]

				except:
					print('----------Warning: Iteration ' + str(i) + '&' + str(j) + ' has initialization error!----------')
					COP[i][j] = None
					Q_c[i][j] = None	

		return COP, Q_c



# =============================================================================
# Load Desiccant Cooling Model and Run Parametric Study Cases
# =============================================================================
class runDesiccantCoolingModel():

	def __init__(self, T_out, RH_out, T_in, RH_in, latentCooling, modelPath, theta_o, V_a, COP_h, r, P_dw, P_fan,
		deltaT, deltaw, phi_min, deltaQ, deltaWc, W_fan, epsilon, epsilon_b, t_g, n, H, d_p, r_R, gCHXisOn):
		"""
		Inputs:
		# 	modelPath: str, add .fmu file path to the python environment
		#   theta_o: float, dimensionless temperature for outlet process air in the DW
		#   V_a: float, volumetric flow rate of moist air, unit: m3/s
		#   COP_h: float, COP of the electric heater
		#	r: float, ratio of mass flow rate of exhaust air in ERW to supply air
		#	deltaT: float, temperature threshold control on/off for the GCHX system, unit: K
		#	deltaw: float, humidity threshold control on/off for the desiccant wheel
		#	phi_min: float, minimum absolute humidity for supply air
		#	deltaQ = float, threshold control on total cooling load, unit: W
		#	deltaWc = float, threshold control on minimum cooling work, unit: W
		#	W_fan = float, fan power, unit: W
		#	epsilon: float, performance of HX in sensible heat recovery
		#	epsilon_b: float, performance of ground heat exchanger in sensible heat recovery
		#	t_g: float, annual average ground temperature, unit: K
		#	n: float, number of parallel boreholes
		#	H: float, depth of the borehole, unit: m
		#	d_p: float, diameter of pipe, unit: m
		#	r_R: float, ratio of thermal resistance of borehole to the ground
		#	COP_c: float, COP of cooling coil on substitution of GHX
		"""
		self.T_out = T_out
		self.RH_out = RH_out
		self.T_in = T_in
		self.RH_in = RH_in
		self.q_lat = latentCooling
		self.model = modelPath
		self.theta_o = theta_o
		self.V_a = V_a
		self.COP_h = COP_h
		self.r = r
		self.deltaT = deltaT
		self.deltaw = deltaw
		self.phi_min = phi_min
		self.deltaQ = deltaQ
		self.deltaWc = deltaWc
		self.W_fan = W_fan
		self.epsilon = epsilon
		self.epsilon_b = epsilon_b
		self.t_g = t_g
		self.n = n
		self.H = H
		self.d_p = d_p
		self.r_R = r_R
		self.P_dw = P_dw
		self.P_fan = P_fan
		self.gCHXisOn = gCHXisOn

	def runModelicaModel(self): 

		# Prepare the result output
		COP = np.zeros((len(self.T_out), len(self.RH_out)))
		Q_c = np.zeros((len(self.T_out), len(self.RH_out)))
		T_lc = np.zeros((len(self.T_out), len(self.RH_out)))
		T_ec = np.zeros((len(self.T_out), len(self.RH_out)))

		# Load the model
		model = load_fmu(self.model)

		# Run quasi-steady model for each time

		for i in range(len(self.T_out)):
			for j in range(len(self.RH_out)):
				# Prepare the inputs
				var = ['t_o', 'rh_o', 't_i', 'rh_i', 'q_lat']
				t = [0.0]
				t_o = [self.T_out[i]]
				rh_o = [self.RH_out[j]]
				t_i = [self.T_in[i]]
				rh_i = [self.RH_in[j]]
				q_lat = [self.q_lat[i][j]]

				# Reset model inputs and parameters
				model.reset()

				# Model parameter setting
				theta_o = self.theta_o
				model.set('theta_o', theta_o)
				V_a = self.V_a
				model.set('V_a', V_a)
				COP_h = self.COP_h
				model.set('COP_h', COP_h)
				r = self.r
				model.set('r', r)
				deltaT = self.deltaT
				model.set('deltaT', deltaT)
				deltaw = self.deltaw
				model.set('deltaw', deltaw)
				phi_min = self.phi_min
				model.set('phi_min', phi_min)
				deltaQ = self.deltaQ
				model.set('deltaQ', deltaQ)
				deltaWc = self.deltaWc
				model.set('deltaWc', deltaWc)
				W_fan = self.W_fan
				model.set('W_fan', W_fan)
				epsilon = self.epsilon
				model.set('epsilon', epsilon)
				epsilon_b = self.epsilon_b
				model.set('epsilon_b', epsilon_b)
				gCHXisOn = self.gCHXisOn
				model.set('gCHXisOn', gCHXisOn)
				t_g = self.t_g
				model.set('t_g', t_g)
				n = self.n
				model.set('n', n)
				H = self.H
				model.set('H', H)
				d_p = self.d_p
				model.set('d_p', d_p)
				r_R = self.r_R
				model.set('r_R', r_R)
				P_fan = self.P_fan
				model.set('P_fan', P_fan)
				P_dw = self.P_dw
				model.set('P_dw', P_dw)


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
					COP[i][j] = res['COP'][-1]
					Q_c[i][j] = res['q_c'][-1]

				except:
					print('----------Warning: Iteration ' + str(i) + '&' + str(j) + ' has initialization error!----------')
					COP[i][j] = None
					Q_c[i][j] = None
					
		return COP, Q_c



# =============================================================================
# Load Membrane Cooling Model and Run Parametric Study Cases
# =============================================================================

class runMembraneCoolingModel():

	def __init__(self, T_out, RH_out, T_in, RH_in, latentCooling, modelPath, epsilon_s, epsilon_l, V_a, eta_is, 
		P_stage, r, deltaT, deltaw, deltaP, phi_min, deltaQ, deltaWc, W_fan, epsilon, epsilon_b, P_fan,
		t_g, n, H, d_p, r_R, A_m, n_m, h, p_w, gCHXisOn):
		"""
		Inputs:
		# 	csvFilePath: str, add .csv file path to the python environment
		# 	modelPath: str, add .fmu file path to the python environment
		# 	idx_start: int, index of start time for simulation
		# 	idx_end: int, index of end time for simulation
		#   epsilon_s: float, performance of ERW in sensible heat recovery
		#   epsilon_l: float, performance of ERW in latent heat recovery
		#   V_a: float, volumetric flow rate of moist air, unit: m3/s
		#   eta_is: float, performance of vacuum pump
		#   P_sweep: float, vacuum pressure at the sweep side of the membrane, Unit: Pa, 611.7 <= P_sweep <= 10M
		#   P_stage: float, stage pressure of the pump, unit: Pa
		#   COP_cond: float, COP of the condenser
		#	r: float, ratio of mass flow rate of exhaust air in ERW to supply air
		#	deltaT: float, temperature threshold control on/off for the ERV system, unit: K
		#	deltaw: float, humidity threshold control on/off for the ERV system
		#	deltaP: float, pressure threshold control on/off for the ERV system, unit: Pa
		#	phi_min: float, minimum absolute humidity for supply air
		#	deltaQ = float, threshold control on total cooling load, unit: W
		#	deltaWc = float, threshold control on minimum cooling work, unit: W
		#	W_fan = float, fan power, unit: W
		#	epsilon: float, performance of GCHX in sensible heat recovery
		#	epsilon_b: float, performance of ground heat exchanger in sensible heat recovery
		#	t_g: float, annual average ground temperature, unit: K
		#	n: float, number of parallel boreholes
		#	H: float, depth of the borehole, unit: m
		#	d_p: float, diameter of pipe, unit: m
		#	r_R: float, ratio of thermal resistance of borehole to the ground
		"""
		self.T_out = T_out
		self.RH_out = RH_out
		self.T_in = T_in
		self.RH_in = RH_in
		self.q_lat = latentCooling
		self.model = modelPath
		self.epsilon_s = epsilon_s
		self.epsilon_l = epsilon_l
		self.V_a = V_a
		self.eta_is = eta_is
		self.P_stage = P_stage
		self.r = r
		self.deltaT = deltaT
		self.deltaw = deltaw
		self.deltaP = deltaP
		self.phi_min = phi_min
		self.deltaQ = deltaQ
		self.deltaWc = deltaWc
		self.W_fan = W_fan
		self.epsilon = epsilon
		self.epsilon_b = epsilon_b
		self.t_g = t_g
		self.n = n
		self.H = H
		self.d_p = d_p
		self.r_R = r_R
		self.A_m = A_m
		self.n_m = n_m
		self.h = h
		self.p_w = p_w
		self.P_fan = P_fan
		self.gCHXisOn = gCHXisOn


	def runModelicaModel(self): 

		# Prepare the result output
		COP = np.zeros((len(self.T_out), len(self.RH_out)))
		Q_c = np.zeros((len(self.T_out), len(self.RH_out)))
		T_supply = np.zeros((len(self.T_out), len(self.RH_out)))
		W_tot = np.zeros((len(self.T_out), len(self.RH_out)))

		# Load the model
		model = load_fmu(self.model)

		# Run quasi-steady model for each time

		for i in range(len(self.T_out)):
			for j in range(len(self.RH_out)):
				# Prepare the inputs
				var = ['t_o', 'rh_o', 't_i', 'rh_i', 'q_lat']
				t = [0.0]
				t_o = [self.T_out[i]]
				rh_o = [self.RH_out[j]]
				t_i = [self.T_in[i]]
				rh_i = [self.RH_in[j]]
				q_lat = [self.q_lat[i][j]]

				# Reset model inputs and parameters
				model.reset()

				# Model parameter setting
				epsilon_s = self.epsilon_s
				model.set('epsilon_s', epsilon_s)
				epsilon_l = self.epsilon_l
				model.set('epsilon_l', epsilon_l)
				V_a = self.V_a
				model.set('V_a', V_a)
				eta_is = self.eta_is
				model.set('eta_is', eta_is)
				P_stage = self.P_stage
				model.set('P_stage', P_stage)
				r = self.r
				model.set('r', r)
				deltaT = self.deltaT
				model.set('deltaT', deltaT)
				deltaw = self.deltaw
				model.set('deltaw', deltaw)
				deltaP = self.deltaP
				model.set('deltaP', deltaP)
				phi_min = self.phi_min
				model.set('phi_min', phi_min)
				deltaQ = self.deltaQ
				model.set('deltaQ', deltaQ)
				deltaWc = self.deltaWc
				model.set('deltaWc', deltaWc)
				W_fan = self.W_fan
				model.set('W_fan', W_fan)
				epsilon = self.epsilon
				model.set('epsilon', epsilon)
				epsilon_b = self.epsilon_b
				model.set('epsilon_b', epsilon_b)
				t_g = self.t_g
				model.set('t_g', t_g)
				n = self.n
				model.set('n', n)
				H = self.H
				model.set('H', H)
				d_p = self.d_p
				model.set('d_p', d_p)
				r_R = self.r_R
				model.set('r_R', r_R)
				A_m = self.A_m
				model.set('A_m', A_m)
				n_m = self.n_m
				model.set('n_m', n_m)
				h = self.h
				model.set('h', h)
				p_w = self.p_w
				model.set('p_w', p_w)
				P_fan = self.P_fan
				model.set('P_fan', P_fan)
				gCHXisOn = self.gCHXisOn
				model.set('gCHXisOn', gCHXisOn)

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
					COP[i][j] = res['COP'][-1]
					Q_c[i][j] = res['q_c'][-1]

				except:
					print('----------Warning: Iteration ' + str(i) + '&' + str(j) + ' has initialization error!----------')
					COP[i][j] = None
					Q_c[i][j] = None

		return COP, Q_c



# =============================================================================
# Load Membrane Cooling With Membrane Model and Run Parametric Study Cases
# =============================================================================

class runMembraneCoolingWithMembModel():

	def __init__(self, T_out, RH_out, T_in, RH_in, latentCooling, modelPath, V_a, epsilon, r, deltaw, deltaP,
		phi_min, deltaQ, deltaWc, W_fan, A_m, n_m, h, p_w, epsilon_b, t_g, n, H, d_p, 
		r_R, COP_h, deltaT, eta_is, P_fan, gCHXisOn):
		"""
		Inputs:
		# 	csvFilePath: str, add .csv file path to the python environment
		# 	modelPath: str, add .fmu file path to the python environment
		# 	idx_start: int, index of start time for simulation
		# 	idx_end: int, index of end time for simulation
		#   epsilon_s: float, performance of ERW in sensible heat recovery
		#   epsilon_l: float, performance of ERW in latent heat recovery
		#   V_a: float, volumetric flow rate of moist air, unit: m3/s
		#   eta_is: float, performance of vacuum pump
		#   P_sweep: float, vacuum pressure at the sweep side of the membrane, Unit: Pa, 611.7 <= P_sweep <= 10M
		#   P_stage: float, stage pressure of the pump, unit: Pa
		#   COP_cond: float, COP of the condenser
		#	r: float, ratio of mass flow rate of exhaust air in ERW to supply air
		#	deltaT: float, temperature threshold control on/off for the ERV system, unit: K
		#	deltaw: float, humidity threshold control on/off for the ERV system
		#	deltaP: float, pressure threshold control on/off for the ERV system, unit: Pa
		#	phi_min: float, minimum absolute humidity for supply air
		#	deltaQ = float, threshold control on total cooling load, unit: W
		#	deltaWc = float, threshold control on minimum cooling work, unit: W
		#	W_fan = float, fan power, unit: W
		#	epsilon: float, performance of GCHX in sensible heat recovery
		#	epsilon_b: float, performance of ground heat exchanger in sensible heat recovery
		#	t_g: float, annual average ground temperature, unit: K
		#	n: float, number of parallel boreholes
		#	H: float, depth of the borehole, unit: m
		#	d_p: float, diameter of pipe, unit: m
		#	r_R: float, ratio of thermal resistance of borehole to the ground
		"""
		self.T_out = T_out
		self.RH_out = RH_out
		self.T_in = T_in
		self.RH_in = RH_in
		self.q_lat = latentCooling
		self.model = modelPath
		self.V_a = V_a
		self.r = r
		self.deltaw = deltaw
		self.deltaP = deltaP
		self.phi_min = phi_min
		self.deltaQ = deltaQ
		self.deltaWc = deltaWc
		self.W_fan = W_fan
		self.epsilon = epsilon
		self.A_m = A_m
		self.n_m = n_m
		self.h = h
		self.p_w = p_w
		self.epsilon_b = epsilon_b
		self.t_g = t_g
		self.n = n
		self.H = H
		self.d_p = d_p
		self.r_R = r_R
		self.COP_h = COP_h
		self.deltaT = deltaT
		self.eta_is = eta_is
		self.gCHXisOn = gCHXisOn


	def runModelicaModel(self): 

		# Prepare the result output
		COP = np.zeros((len(self.T_out), len(self.RH_out)))
		Q_c = np.zeros((len(self.T_out), len(self.RH_out)))
		T_supply = np.zeros((len(self.T_out), len(self.RH_out)))
		W_tot = np.zeros((len(self.T_out), len(self.RH_out)))

		# Load the model
		model = load_fmu(self.model)

		# Run quasi-steady model for each time

		for i in range(len(self.T_out)):
			for j in range(len(self.RH_out)):
				# Prepare the inputs
				var = ['t_o', 'rh_o', 't_i', 'rh_i', 'q_lat']
				t = [0.0]
				t_o = [self.T_out[i]]
				rh_o = [self.RH_out[j]]
				t_i = [self.T_in[i]]
				rh_i = [self.RH_in[j]]
				q_lat = [self.q_lat[i][j]]

				# Reset model inputs and parameters
				model.reset()

				# Model parameter setting
				V_a = self.V_a
				model.set('V_a', V_a)
				r = self.r
				model.set('r', r)
				deltaw = self.deltaw
				model.set('deltaw', deltaw)
				deltaP = self.deltaP
				model.set('deltaP', deltaP)
				phi_min = self.phi_min
				model.set('phi_min', phi_min)
				deltaQ = self.deltaQ
				model.set('deltaQ', deltaQ)
				deltaWc = self.deltaWc
				model.set('deltaWc', deltaWc)
				W_fan = self.W_fan
				model.set('W_fan', W_fan)
				epsilon = self.epsilon
				model.set('epsilon', epsilon)
				A_m = self.A_m
				model.set('A_m', A_m)
				n_m = self.n_m
				model.set('n_m', n_m)
				h = self.h
				model.set('h', h)
				p_w = self.p_w
				model.set('p_w', p_w)
				epsilon_b = self.epsilon_b
				model.set('epsilon_b', epsilon_b)
				t_g = self.t_g
				model.set('t_g', t_g)
				n = self.n
				model.set('n', n)
				H = self.H
				model.set('H', H)
				d_p = self.d_p
				model.set('d_p', d_p)
				r_R = self.r_R
				model.set('r_R', r_R)
				COP_h = self.COP_h
				model.set('COP_h', COP_h)
				deltaT = self.deltaT
				model.set('deltaT', deltaT)
				eta_is = self.eta_is
				model.set('eta_is', eta_is)
				gCHXisOn = self.gCHXisOn
				model.set('gCHXisOn', gCHXisOn)

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
					COP[i][j] = res['COP'][-1]
					Q_c[i][j] = res['q_c'][-1]

				except:
					print('----------Warning: Iteration ' + str(i) + '&' + str(j) + ' has initialization error!----------')
					COP[i][j] = None
					Q_c[i][j] = None

		return COP, Q_c