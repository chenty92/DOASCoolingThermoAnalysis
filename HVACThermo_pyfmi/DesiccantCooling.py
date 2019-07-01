"""
@author: Tianyi Chen
Version 1.0, 09/06/2018, Framework completed.
"""

from pyfmi import load_fmu
import numpy as np
import EPlusOutput as epo

# =============================================================================
# Load Modelica Model and Run Parametric Study Cases
# =============================================================================

class runDesiccantCoolingModel():

	def __init__(self, csvFilePath, modelPath, idx_start, idx_end, theta_o, V_a, COP_h, r, P_dw, P_fan,
		deltaT, deltaw, phi_min, deltaQ, deltaWc, W_fan, epsilon, epsilon_b, t_g, n, H, d_p, r_R, gCHXisOn = 1):
		"""
		Inputs:
		# 	csvFilePath: str, add .csv file path to the python environment
		# 	modelPath: str, add .fmu file path to the python environment
		# 	idx_start: int, index of start time for simulation
		# 	idx_end: int, index of end time for simulation
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
		self.path = csvFilePath
		self.model = modelPath
		self.start = idx_start
		self.end = idx_end
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
		"""
		Outputs:
		#   COP: list[float], time series solutions of COP for the desiccant system in kW, len(COP) = idx_end - idx_start
		#   W_tot: list[float], time series solutions of total work for the desiccant system in kW, len(W_tot) = idx_end - idx_start
		#   Q_c: list[float], time series solutions of cooling load for the desiccant system in kW, len(Q_c) = idx_end - idx_start
		"""

		# Read the EPlus output file
		data = epo.EPlusOutput(self.path)
		time, T_out, RH_out, T_in, RH_in, cooling, sensibleCooling, latentCooling = data.readFile()

		# Prepare the result output
		COP = [0] * (self.end - self.start)
		T_supply = [0] * (self.end - self.start)
		W_tot = [0] * (self.end - self.start)
		Q_c = [0] * (self.end - self.start)


		# Load the model
		model = load_fmu(self.model, log_level = 7)

		# Run quasi-steady model for each time

		for i in range(self.start, self.end):
			# Prepare the inputs
			var = ['t_o', 'rh_o', 't_i', 'rh_i', 'q_lat']
			t = [0.0]
			t_o = [T_out[i]]
			rh_o = [RH_out[i]]
			t_i = [T_in[i]]
			rh_i = [RH_in[i]]
			q_lat = [latentCooling[i]]

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
				res = model.simulate(start_time = 0, final_time = 1.0, input = input_object, options = simulate_opts)

				# Results
				COP[i - self.start] = res['COP'][-1]
				Q_c[i - self.start] = res['q_c'][-1]

			except:
				print('----------Warning: Iteration ' + str(i) + ' has initialization error!----------')
				COP[i - self.start] = None
				Q_c[i - self.start] = None

		return COP, Q_c

			# print i
			# print "DW inlet temperature =", res['desiccantWheel.t_1'][-1]
			# print "DW inlet humidity ratio =", res['desiccantWheel.w_1'][-1]
			# print "DW outlet humidity ratio =", res['desiccantWheel.w_2'][-1]
			# print "DW dehumidified air temperature =", res['desiccantWheel.t_2'][-1]
			# print "DW regeneration air temperature =", res['desiccantWheel.T_3'][-1]
			# print "DW heating work =", res['desiccantWheel.w_h'][-1]
			# print "GCHX inlet temperature =", res['gCHX.t_1'][-1]
			# print "GCHX pump work =", res['gCHX.w_c'][-1]
