"""
@author: Tianyi Chen
Version 1.0, 12/17/2018, Framework completed.
"""

from pyfmi import load_fmu
import numpy as np
import EPlusOutput as epo

# =============================================================================
# Load Modelica Model and Run Parametric Study Cases
# =============================================================================

class runMembraneCoolingWithMembModel():

	def __init__(self, csvFilePath, modelPath, idx_start, idx_end, V_a, epsilon, r, deltaw, deltaP,
		phi_min, deltaQ, deltaWc, W_fan, A_m, n_m, h, p_w, epsilon_b, t_g, n, H, d_p, 
		r_R, COP_h, deltaT, eta_is, P_fan, gCHXisOn = 1):

		self.path = csvFilePath
		self.model = modelPath
		self.start = idx_start
		self.end = idx_end
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
				res = model.simulate(start_time = 0, final_time = 1.0, input = input_object, options = simulate_opts)

				# Results
				COP[i - self.start] = res['COP'][-1]
				Q_c[i - self.start] = res['q_c'][-1]

				# print res['membraneWithMembrane.T_4'][-1]
				# print res['membraneWithMembrane.T_5'][-1]
				# print res['membraneWithMembrane.T_5min'][-1]
				# print res['membraneWithMembrane.P_sweep'][-1]
				# print res['membraneWithMembrane.P_stage'][-1]
				# print res['membraneWithMembrane.W_h'][-1]
				# print res['membraneWithMembrane.W_stage'][-1]


			except:
				print('----------Warning: Iteration ' + str(i) + ' has initialization error!----------')
				COP[i - self.start] = None
				Q_c[i - self.start] = None	

		return COP, Q_c
