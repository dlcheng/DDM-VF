# This routine returns the transfer function 

import vals
import numpy as np
from scipy import interpolate

def prim_tk(k):
	x = 1
	if vals.Ps_flag == 1:
		x = 1
	if vals.Ps_flag == 2:
		x = tk_bbks(k)
	if vals.Ps_flag == 3:
		pass                     # return E&H transfer function
	if vals.Ps_flag == 4:
		x = tk_file(k)           # return tk from file
	return x

def tk_bbks(k):
	x = 1.0
	T_sig = vals.Omega_m * vals.h0
	T_sig *= np.exp(-1.0 * vals.Omega_b - np.sqrt(2.0 * vals.h0) * vals.Omega_b / vals.Omega_m)
	q = k / T_sig
	a0 = 2.34
	a1 = 3.89
	a2 = 16.19
	a3 = 5.46
	a4 = 6.71
	if q >= 1e-8:                
		b1 = 1.0 
		b1 += a1 * q
		b1 += a2 * a2 * q * q
		b1 += a3 * a3 * a3 * q * q * q
		b1 += a4 * a4 * a4 * a4 * q * q * q * q
		b1 = 1.0 / np.power(b1, 0.25)
		b2 = np.log(1.0 + a0 * q) / a0 / q
		x = b1 * b2
	return x

def init_tk_from_file():
	file_name = './TK/transf0.dat'
	vals.Tb_tk = np.loadtxt(file_name)
	fb = vals.Omega_b / vals.Omega_m
#	fb = 0.166                   # for the WMAP shape of Jing's simulations
	vals.Tb_tk[:,1] = (1-fb) * vals.Tb_tk[:,1] + fb * vals.Tb_tk[:,2]
	vals.Sp_tk = interpolate.InterpolatedUnivariateSpline(np.log10(vals.Tb_tk[:,0]), vals.Tb_tk[:,1]/vals.Tb_tk[0,1])

def tk_file(k):    
    if k >= vals.Tb_tk[0,0] and k<= vals.Tb_tk[-1,0]:
    	return vals.Sp_tk(np.log10(k))
    if k < vals.Tb_tk[0,0]:
    	return 1.0
    if k > vals.Tb_tk[-1,0]:
    	return 0.0

