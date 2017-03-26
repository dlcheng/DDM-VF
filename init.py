# This routine does the initialization of cross-module variables, power spectrum
# as well as the spline function of sigma2 v.s. M with both of them in log10 scale
import vals
import growth
import sigma2
import constants
import tk
import mf

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

def set_cosmological_parameters(Omega_m, Omega_v, Omega_b, h0):
	vals.Omega_m = Omega_m
	vals.Omega_v = Omega_v
	vals.Omega_b = Omega_b
	vals.h0 = h0

def set_power_parameters(ns, sigma8, Ps_flag):
	vals.ns = ns
	vals.sigma8 = sigma8
	vals.Ps_flag = Ps_flag

def init_M_for_MF(Min, Max, Num):
	vals.Tb_M = np.linspace(np.log10(Min), np.log10(Max), Num)
	vals.Tb_M = np.power(10, vals.Tb_M)

def init_growth_factor():
	vals.Agr = 1.0 / growth.unnorm_growth(0.0)	

def init_power_spectrum():
	if vals.Ps_flag == 4:
		tk.init_tk_from_file()
	sigma2.norm_ps()

def init_logsigma2_slpine():
	Num = vals.Tb_M.shape[0]
	vals.Tb_sigma2 = np.zeros((Num,))               # the first step is to prepare the data
	for i in range(Num):
		vals.Tb_sigma2[i] = sigma2.calculate_sigma2(vals.Tb_M[i])
	vals.Sp_logsigma2_logm = interpolate.InterpolatedUnivariateSpline(np.log10(vals.Tb_M), np.log10(vals.Tb_sigma2))

def init_boxsize(box):
	vals.boxsize = box

def init_halo_mf_table(Min, Max, Num):
	vals.Tb_MF = np.zeros((Num, 4))
	dis = np.log10(Max/Min)/Num
	vals.Tb_MF[:,0] = np.linspace(np.log10(Min), np.log10(Max)-dis, Num)
	vals.Tb_MF[:,1] = np.linspace(np.log10(Min)+dis, np.log10(Max), Num)
	vals.Tb_MF[:,2] = (vals.Tb_MF[:,0]+vals.Tb_MF[:,1])/2.0
	vals.Tb_MF[:,3] = mf.mass_function(np.power(10, vals.Tb_MF[:,2])) * dis * np.power(vals.boxsize, 3)

def init_ddm_param(vk, fd):
	vals.vk = vk
	vals.fd = fd
