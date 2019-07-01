"""
@author: Tianyi Chen
Version 1.0, 05/02/2018, Framework completed.
Version 1.1, 05/10/2018, Parameter sets expanded with control threshold values.
Version 1.2, 05/15/2018, Baseline Chiller model with COP fitting curve in Modelica.
Version 1.3, 09/06/2018, Membrane and desiccant cooling models are ready.
"""

from HVACThermo_pyfmi import EPlusOutput as epo
from HVACThermo_pyfmi import BaselineChiller as bc
from HVACThermo_pyfmi import MembraneCooling as mc
from HVACThermo_pyfmi import DesiccantCooling as dc
from HVACThermo_pyfmi import ERWChiller as ec
from HVACThermo_pyfmi import HRWChiller as hc
from HVACThermo_pyfmi import MembraneCoolingWithMemb as mcm
# from HVACThermo_pyfmi import MembraneCooling_SG as mc_sg
import numpy as np
import matplotlib.pyplot as plt
import csv


# =============================================================================
# General Settings
# =============================================================================
##### Read weather/load data file
# # Case 1: BOS Low-energy Building Prototype
# root = 'C:/temp/LowE_Off/'
# name = 'loE'
# Case 2: Passive Strategies
root = 'D:/tianyich/Dropbox (MIT)/PhD Work/PassiveStrategies4Cooling/'
name = 'Boston_Hybrid_HighTM'
# Combine file path
csvFilePath = root + name + '.csv'

# # Month: July
# idx_start = 4344
# idx_end = 5087
# # Week: 07/10 - 07/17
# idx_start = 4560
# idx_end = 4727
# # Day of 08/15
# idx_start = 5424
# idx_end = 5447
# Whole Year
idx_start = 0
idx_end = 8760
# # Test
# idx_start = 1
# idx_end = 2

##### Parameter Settings
# Ventilation rate
TFA = 1750						# Unit: m2
V_a = 0.55 * 1e-3 * TFA 		# Unit: m3/s

# TFA = 100						# Unit: m2
# V_a = 0.55 * 1e-3 * TFA 		# Unit: m3/s

# Supply/exhaust fan
W_fan = 1.5 * V_a * 1e3			# Unit: W

# HRW / HX
epsilon = 0.64
P_HRW = 241

# HRW / HX: dcs
epsilon_dcs = 0.80
P_HRW_dcs = 85

# ERW
epsilon_s = 0.64
epsilon_l = 0.58
P_ERW = 241

# ERW: mcs
epsilon_s_mcs = 0.80
epsilon_l_mcs = 0.77
P_ERW_mcs = 85

# Chiller
unit = 1.0
V_w = 0.12113 / unit 			# Unit: m3/s
t_lcmin = 4.44					# Unit: C

# Desiccant
theta_o = 1.75
COP_h = 4.8
P_dw = 120

# Membrane
eta_is = 0.80
P_stage = 2000
A_m_mccs = 5.0
A_m_mcms = 5.0
n_m = 500.0
h = 0.005
p_w = 6.2e-7
COP_h = 4.8

# GHE
epsilon_b = 0.60
t_g = 10.0
n = 1.0
H = 50.0
d_p = 0.025
r_R = 0.42
############################################################################


# =============================================================================
# Case Study for Chiller Model  
# =============================================================================
"""
Input paraters
# csvFilePath: str, add .csv file path to the python environment
# modelPath: str, add .fmu file path to the python environment
# idx_start: int, index of start time for simulation
# idx_end: int, index of end time for simulation
# V_a: float, volumetric flow rate of moist air, unit: m3/s
# V_w: float, volumetric flow rate of chilled water, unit: m3/s, V_w = refVw / floor
# epsilon: float, performance of HX in chiller
# unit: float, number of number of parallel terminals
# phi_min: float, minimum absolute humidity for supply air
# t_lcmin: float, minimum leaving chilled water temperature, unit: C
# deltaCOP: float, threshold control on chiller
# deltaQ = float, threshold control on total cooling load, unit: W
# deltaWc = float, threshold control on minimum cooling work, unit: W
# W_fan = float, fan power, unit: W
"""
########################################################################################################
# Modelica Model Input Parameters
modelPath = '../Working Directory/HVACThermo_BaselineChiller.fmu'
########################################################################################################

# Run Simulation for Chiller
model = bc.runBaselineChillerModel(csvFilePath = csvFilePath, modelPath = modelPath, idx_start = idx_start, idx_end = idx_end, 
	V_a = V_a, V_w = V_w, epsilon = epsilon, phi_min = 0.1, t_lcmin = t_lcmin, deltaCOP = 1e-2, deltaQ = 1.0, 
	deltaWc = 1.0, W_fan = W_fan)
COP_bcs, Q_c_bcs = model.runModelicaModel()



# =============================================================================
# Case Study for Desiccant Cooling Model
# =============================================================================
"""
Input paraters
# csvFilePath: str, add .csv file path to the python environment
# modelPath: str, add .fmu file path to the python environment
# idx_start: int, index of start time for simulation
# idx_end: int, index of end time for simulation
# theta_o: float, dimensionless temperature for outlet process air in the DW
# COP_h: float, COP of the electric heater
# r: float, ratio of mass flow rate of exhaust air in ERW to supply air
# deltaT: float, temperature threshold control on/off for the GCHX system, unit: K
# deltaw: float, humidity threshold control on/off for the desiccant wheel
# phi_min: float, minimum absolute humidity for supply air
# deltaQ = float, threshold control on total cooling load, unit: W
# deltaWc = float, threshold control on minimum cooling work, unit: W
# W_fan = float, fan power, unit: W
# epsilon: float, performance of HX in sensible heat recovery
# epsilon_b: float, performance of ground heat exchanger in sensible heat recovery
# t_g: float, annual average ground temperature, unit: K
# n: float, number of parallel boreholes
# H: float, depth of the borehole, unit: m
# r_p: float, diameter of pipe, unit: m
# COP_c: float, COP of cooling coil on substitution of GHX
"""
#######################################################################################################
# Modelica Model Input Parameters
modelPath = '../Working Directory/HVACThermo_DesiccantCooling.fmu'
########################################################################################################

# Run Simulation
model = dc.runDesiccantCoolingModel(csvFilePath = csvFilePath, modelPath = modelPath, idx_start = idx_start, idx_end = idx_end, 
		theta_o = theta_o, V_a = V_a, COP_h = COP_h, r = 1.0, deltaT = 1, deltaw = 1e-5, phi_min = 0.1, 
		deltaQ = 1.0, deltaWc = 1.0, W_fan = W_fan, epsilon = epsilon_dcs, epsilon_b = epsilon_b, t_g = t_g, 
		n = n, H = H, d_p = d_p, r_R = r_R, P_dw = P_dw, P_fan = P_HRW_dcs)
COP_dcs, Q_c_dcs = model.runModelicaModel()



# =============================================================================
# Case Study for Membrane Cooling Model
# =============================================================================
"""
Input paraters
# csvFilePath: str, add .csv file path to the python environment
# modelPath: str, add .fmu file path to the python environment
# idx_start: int, index of start time for simulation
# idx_end: int, index of end time for simulation
# epsilon_s: float, performance of ERW in sensible heat recovery
# epsilon_l: float, performance of ERW in latent heat recovery
# V_a: float, volumetric flow rate of moist air, unit: m3/s
# eta_is: float, performance of vacuum pump
# P_sweep: float, vacuum pressure at the sweep side of the membrane, Unit: Pa, 611.7 <= P_sweep <= 10M
# P_stage: float, stage pressure of the pump, unit: Pa
# COP_cond: float, COP of the condenser
# epsilon: float, performance of GCHX in sensible heat recovery
# epsilon_b: float, performance of ground heat exchanger in sensible heat recovery
# t_g: float, annual average ground temperature, unit: K
# n: float, number of parallel boreholes
# H: float, depth of the borehole, unit: m
# d_p: float, diameter of pipe, unit: m
# r_R: float, ratio of thermal resistance of borehole to the ground
"""
########################################################################################################
# Modelica Model Input Parameters
modelPath = '../Working Directory/HVACThermo_MembraneCooling.fmu'
A_m = A_m_mccs
########################################################################################################

# Run Simulation
model = mc.runMembraneCoolingModel(csvFilePath = csvFilePath, modelPath = modelPath, idx_start = idx_start, idx_end = idx_end, 
		epsilon_s = epsilon_s_mcs, epsilon_l = epsilon_l_mcs, V_a = V_a, eta_is = eta_is, P_stage = P_stage, P_fan = P_ERW_mcs,
		r = 1.0, deltaT = 1, deltaw = 1e-5, deltaP = 1.0, phi_min = 0.1, deltaQ = 1.0, deltaWc = 1.0, 
		W_fan = W_fan, epsilon = epsilon_s_mcs, epsilon_b = epsilon_b, t_g = t_g, n = n, H = H, d_p = d_p, r_R = r_R,
		A_m = A_m, n_m = n_m, h = h, p_w = p_w)
COP_mcs, Q_c_mcs = model.runModelicaModel()



# =============================================================================
# Case Study for Membrane Cooling With Membrane Model  
# =============================================================================
"""
Input paraters
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
########################################################################################################
# Modelica Model Input Parameters
modelPath = '../Working Directory/HVACThermo_MembraneCoolingWithMemb.fmu'
A_m = A_m_mcms
########################################################################################################

# Run Simulation for HRWChiller
model = mcm.runMembraneCoolingWithMembModel(csvFilePath = csvFilePath, modelPath = modelPath, idx_start = idx_start, idx_end = idx_end, 
	V_a = V_a, epsilon = epsilon, phi_min = 0.1, deltaQ = 1.0, deltaWc = 1.0, r = 1.0, eta_is = eta_is, P_fan = P_HRW,
	W_fan = W_fan, A_m = A_m, n_m = n_m, h = h, p_w = p_w, deltaw = 1e-5, deltaP = 1.0, epsilon_b = epsilon_b, t_g  = t_g, n = n, 
	H = H, d_p = d_p, r_R = r_R, COP_h = COP_h, deltaT = 1.0)
COP_mcms, Q_c_mcms = model.runModelicaModel()



# =============================================================================
# Case Study for ERWChiller Model  
# =============================================================================
"""
Input paraters
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
########################################################################################################
# Modelica Model Input Parameters
modelPath = '../Working Directory/HVACThermo_ERWChiller.fmu'
########################################################################################################

# Run Simulation for ERWChiller
model = ec.runERWChillerModel(csvFilePath = csvFilePath, modelPath = modelPath, idx_start = idx_start, idx_end = idx_end, 
	V_a = V_a, V_w = V_w, epsilon = epsilon, epsilon_s = epsilon_s, epsilon_l = epsilon_l, phi_min = 0.1, t_lcmin = t_lcmin, 
	deltaCOP = 1e-2, deltaQ = 1.0, deltaWc = 1.0, r = 1.0, deltaT = 1, deltaw = 1e-5, deltaP = 1.0, W_fan = W_fan, P_fan = P_HRW)
COP_ecs, Q_c_ecs = model.runModelicaModel()



# =============================================================================
# Case Study for HRWChiller Model  
# =============================================================================
"""
Input paraters
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
########################################################################################################
# Modelica Model Input Parameters
modelPath = '../Working Directory/HVACThermo_HRWChiller.fmu'
########################################################################################################

# Run Simulation for HRWChiller
model = hc.runHRWChillerModel(csvFilePath = csvFilePath, modelPath = modelPath, idx_start = idx_start, idx_end = idx_end, 
	V_a = V_a, V_w = V_w, epsilon = epsilon, phi_min = 0.1, t_lcmin = t_lcmin, P_fan = P_HRW,
	deltaCOP = 1e-2, deltaQ = 1.0, deltaWc = 1.0, r = 1.0, deltaT = 1, W_fan = W_fan)
COP_hcs, Q_c_hcs = model.runModelicaModel()


####################################################
# # Write into a .csv file
data = epo.EPlusOutput(csvFilePath = csvFilePath)
time, T_out, RH_out, T_in, RH_in, cooling, sensibleCooling, latentCooling = data.readFile()
cols = ['Date/Time', 'T_out (C)', 'RH_out', 'T_in (C)', 'RH_in', 'Q_lat (kW)', 'Q_c (kW)', 'COP_bcs', 'COP_dcs', 'COP_mcs', 'COP_mcms', 'COP_ecs', 'COP_hcs']
resultFilePath = root + name + '_res.csv'
with open(resultFilePath, 'wb') as f:
	writer = csv.writer(f)
	writer.writerow(cols)
	for i in range(len(COP_bcs)):
		if Q_c_bcs[i] <= 1e-3 or Q_c_mcs[i] <= 1e-3 or Q_c_dcs[i] <= 1e-3 or Q_c_mcms[i] <= 1e-3:
			row = [time[idx_start+i], T_out[idx_start+i], RH_out[idx_start+i], T_in[idx_start+i], RH_in[idx_start+i], latentCooling[idx_start + i]/3.6e6] + [0.0] * (len(cols) - 6)
		else:
			row = [time[idx_start+i], T_out[idx_start+i], RH_out[idx_start+i], T_in[idx_start+i], RH_in[idx_start+i], latentCooling[idx_start + i]/3.6e6, Q_c_bcs[i], COP_bcs[i], COP_dcs[i], COP_mcs[i], COP_mcms[i], COP_ecs[i], COP_hcs[i]]
		writer.writerow(row)
f.close()









# =============================================================================
# Result Analysis
# =============================================================================

# print "COP of baseline chiller =", COP_bcs
# print "COP of membrane cooling", COP_mcs
# print "COP of membrane cooling with Membrane", COP_mcms
# print "COP of desiccant cooling =", COP_dcs

# # Plot COP vs time
# fig = plt.figure()
# plt.rcParams["font.family"] = "Times New Roman"
# time = np.arange(1, len(COP_mcs) + 1)
# plt.plot(time, COP_bcs, 'k-', label = 'Baseline Chiller')
# plt.plot(time, COP_mcs, 'k--', label = 'Membrane Cooling')
# plt.plot(time, COP_dcs, 'k-.', label = 'Desiccant Cooling')
# plt.title('COP of HVAC Cooling Systems in July')
# plt.xlabel('Time (hour)')
# plt.ylabel('COP')
# # plt.ylim(ymax = 15)
# plt.legend()
# # plt.show()
# fig.savefig('Plots/Results/COP_July.png', bbox_inches='tight', dpi=1200)

# # Plot COP vs q_lat
# data = epo.EPlusOutput(csvFilePath = csvFilePath)
# time, T_out, RH_out, T_in, RH_in, cooling, sensibleCooling, latentCooling = data.readFile()
# q_lat = list(np.array(latentCooling[idx_start: idx_end]) / 3.6e3 / 1e3)
# fig = plt.figure()
# plt.rcParams["font.family"] = "Times New Roman"
# plt.plot(Q_c, COP_bcs, 'ko', label = 'ElectricEIRChiller Trane CVHF 2567kW/11.77COP/VSD')
# plt.plot(Q_c, COP_mcs, 'k^', label = 'ERW + Membrane (with condenser)')
# plt.title('COP vs Cooling Load')
# plt.xlabel('Latent Cooling Load (kW)')
# plt.ylabel('COP')
# plt.ylim(ymin = 0)
# plt.legend()
# plt.show()
# fig.savefig('COP vs Latent Cooling Load.png', bbox_inches='tight', dpi=1200)

# # Plot COP vs q_lat
# T_o = T_out[idx_start: idx_end]
# fig = plt.figure()
# plt.rcParams["font.family"] = "Times New Roman"
# plt.plot(T_o, COP_bcs, 'ko', label = 'ElectricEIRChiller Trane CVHF 2567kW/11.77COP/VSD')
# plt.plot(T_o, COP_mcs, 'k^', label = 'ERW + Membrane (with condenser)')
# plt.title('COP vs Outdoor Temperature')
# plt.xlabel('Outdoor Temperature (C)')
# plt.ylabel('COP')
# plt.ylim(ymin = 0)
# plt.legend()
# plt.show()
# fig.savefig('COP vs Outdoor Temperature.png', bbox_inches='tight', dpi=1200)

# # Write into a .csv file
# cols = ['Q_c (kW)', 'COP_bcs', 'COP_dcs', 'COP_mcs', 'COP_mcm', 'COP_ecs', 'COP_hcs']
# resultFilePath = root + 'results/result_' + name + '.csv'
# with open(resultFilePath, 'wb') as f:
# 	writer = csv.writer(f)
# 	writer.writerow(cols)
# 	for i in range(len(COP_bcs)):
# 		if Q_c_bcs[i] <= 1e-3 or Q_c_dcs[i] <= 1e-3 or Q_c_mcs[i] <= 1e-3 or Q_c_mcm[i] <= 1e-3:
# 			row = [0.0] * len(cols)
# 		else:
# 			row = [Q_c_bcs[i], COP_bcs[i], COP_dcs[i], COP_mcs[i], COP_mcm[i], COP_ecs[i], COP_hcs[i]]
# 		writer.writerow(row)
# f.close()


# # Write into a .csv file
# cols = ['Q_c (kW)', 'COP_bcs', 'COP_dcs', 'COP_mcs', 'COP_mcm']
# resultFilePath = 'result_' + name + '.csv'
# with open(resultFilePath, 'wb') as f:
# 	writer = csv.writer(f)
# 	writer.writerow(cols)
# 	for i in range(len(COP_bcs)):
# 		if Q_c_bcs[i] <= 1e-3 or Q_c_mcs[i] <= 1e-3 or Q_c_mcm[i] <= 1e-3 or Q_c_dcs[i] <= 1e-3:
# 			row = [0.0] * len(cols)
# 		else:
# 			row = [Q_c_bcs[i], COP_bcs[i], COP_dcs[i], COP_mcs[i], COP_mcm[i]]
# 		writer.writerow(row)
# f.close()

# # Write into a .csv file
# cols = ['Q_c (kW)', 'COP_bcs', 'COP_dcs', 'COP_mcs']
# resultFilePath = 'result_' + name + '.csv'
# with open(resultFilePath, 'wb') as f:
# 	writer = csv.writer(f)
# 	writer.writerow(cols)
# 	for i in range(len(COP_bcs)):
# 		if Q_c_bcs[i] <= 1e-3 or Q_c_mcs[i] <= 1e-3 or Q_c_dcs[i] <= 1e-3:
# 			row = [0.0] * len(cols)
# 		else:
# 			row = [Q_c_bcs[i], COP_bcs[i], COP_dcs[i], COP_mcs[i]]
# 		writer.writerow(row)
# f.close()


#####################################################

# # =============================================================================
# # Case Study for Membrane Cooling With Sweeping Gas Model  
# # =============================================================================
# """
# Input paraters
# # 	csvFilePath: str, add .csv file path to the python environment
# # 	modelPath: str, add .fmu file path to the python environment
# # 	idx_start: int, index of start time for simulation
# # 	idx_end: int, index of end time for simulation
# #   V_a: float, volumetric flow rate of moist air, unit: m3/s
# #   V_w: float, volumetric flow rate of water, unit: m3/s
# #   epsilon: float, performance of HX in chiller
# #   epsilon_s: float, performance of ERW in sensible heat recovery
# #   epsilon_l: float, performance of ERW in latent heat recovery
# #	r: float, ratio of mass flow rate of exhaust air in ERW to supply air
# #	unit: float, number of number of parallel terminals
# #	phi_min: float, minimum absolute humidity for supply air
# #	t_lcmin: float, minimum leaving chilled water temperature, unit: C
# #	deltaCOP: float, threshold control on chiller
# #	deltaT: float, temperature threshold control on/off for the ERV system, unit: K
# #	deltaw: float, humidity threshold control on/off for the ERV system
# #	deltaP: float, pressure threshold control on/off for the ERV system, unit: Pa
# #	deltaQ: float, threshold control on total cooling load, unit: W
# #	deltaWc: float, threshold control on minimum cooling work, unit: W
# #	W_fan: float, fan power, unit: W
# """
# ########################################################################################################
# # Modelica Model Input Parameters
# modelPath = '../Working Directory/HVACThermo_MembraneCooling_0SG.fmu'
# ########################################################################################################

# # Run Simulation for HRWChiller
# model = mc_sg.runMembraneCoolingSGModel(csvFilePath = csvFilePath, modelPath = modelPath, idx_start = idx_start, idx_end = idx_end, 
# 	V_a = V_a, V_w = V_w, epsilon = epsilon, phi_min = 0.1, t_lcmin = t_lcmin, deltaCOP = 1e-2, deltaQ = 1.0, deltaWc = 1.0, r = 1.0, 
# 	W_fan = W_fan, A_m = A_m, n_m = n_m, h = h, p_w = p_w, deltaw = 1e-5, deltaP = 1.0, epsilon_b = epsilon_b, t_g  = t_g, n = n, 
# 	H = H, d_p = d_p, r_R = r_R, COP_h = COP_h, deltaT = 1.0, eta_is = eta_is)
# COP_mcs_sg, Q_c_mcs_sg = model.runModelicaModel()
