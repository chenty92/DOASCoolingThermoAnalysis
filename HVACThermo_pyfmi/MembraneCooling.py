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

class runMembraneCoolingModel():

	def __init__(self, csvFilePath, modelPath, idx_start, idx_end, epsilon_s, epsilon_l, V_a, eta_is, 
		P_stage, r, deltaT, deltaw, deltaP, phi_min, deltaQ, deltaWc, W_fan, epsilon, epsilon_b, P_fan,
		t_g, n, H, d_p, r_R, A_m, n_m, h, p_w, gCHXisOn = 1):
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
		self.path = csvFilePath
		self.model = modelPath
		self.start = idx_start
		self.end = idx_end
		self.epsilon_s = epsilon_s
		self.epsilon_l = epsilon_l
		self.V_a = V_a
		self.eta_is = eta_is
		# self.P_sweep = P_sweep
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
		"""
		Outputs:
		#   COP: list[float], time series solutions of COP for the membrane system in kW, len(COP) = idx_end - idx_start
		#   T_supply: list[float], time series solutions of supply air temperature in C, len(T_supply) = idx_end - idx_start
		#   W_tot: list[float], time series solutions of total work for the membrane system in kW, len(W_tot) = idx_end - idx_start
		#   Q_c: list[float], time series solutions of cooling load for the membrane system in kW, len(Q_c) = idx_end - idx_start
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
		model = load_fmu(self.model)

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
			epsilon_s = self.epsilon_s
			model.set('epsilon_s', epsilon_s)
			epsilon_l = self.epsilon_l
			model.set('epsilon_l', epsilon_l)
			V_a = self.V_a
			model.set('V_a', V_a)
			eta_is = self.eta_is
			model.set('eta_is', eta_is)
			# P_sweep = self.P_sweep
			# model.set('P_sweep', P_sweep)
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
				res = model.simulate(start_time = 0, final_time = 1.0, input = input_object, options = simulate_opts)

				# Results
				COP[i - self.start] = res['COP'][-1]
				Q_c[i - self.start] = res['q_c'][-1]

				# print res['membraneWithCondenser.w_c'][-1]
				# print res['enthalpyRecoveryWheel.w_p'][-1]
				# print res['gCHX.w_c'][-1]
				# print res['gCHX.t_1'][-1]
				# print res['gCHX.t_2'][-1]


			except:
				print('----------Warning: Iteration ' + str(i) + ' has initialization error!----------')
				COP[i - self.start] = None
				Q_c[i - self.start] = None			

		return COP, Q_c

			# print res['gCHX.T_3'][-1]
			# print res['gCHX.T_4'][-1]
			# print res['gCHX.m_w'][-1]
			# print res['gCHX.v_w'][-1]
			# print res['gCHX.deltaP'][-1]
			# print res['gCHX.w_c'][-1]
			# print i
			# W_tot[i - self.start] = res['w_tot'][-1]
