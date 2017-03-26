import vals
import numpy as np
import constants
import growth

def mass_function(M):                # return dN/dlog10M
	sigma2 = vals.Sp_logsigma2_logm(np.log10(M))
	sigma2 = np.power(10, sigma2)
	rho_m = constants.rho_crit * vals.Omega_m
	v = 1.686 / np.sqrt(sigma2)
	x = 0.5 * np.log(10) * rho_m / M
	x *= fv_st(v)                   # ST formula
#	x *= fv_tinker(v)
	x *= np.abs(vals.Sp_logsigma2_logm(np.log10(M), nu=1))
	x *= 1.25                       # correction for subhaloes
	x *= mf_ddm_sup(M)              # correction for ddm effect
	return x

def fv_st(v, z=0):
	A = 0.3222
	p = 0.3
	q = 0.707
	v = v / growth.gr(z)
	x = A * np.sqrt(2.0/np.pi)
	x *= np.sqrt(q)*v
	x *= 1.0+np.power(np.sqrt(q)*v, -2.0*p)
	x *= np.exp(-0.5 * q * v * v)
	return x

def fv_tinker(v):
	A = 0.224                       # parameter suggested by Klypin et al.(2014)
	a = 1.67
	b = 1.80
	c = 1.48
	sigma = 1.686 / v
	x = A * (np.power(sigma/b, -a) + 1)
	x *= np.exp(-c/sigma/sigma)
	return x

def mf_ddm_sup(M):
	M_f = 4.0*np.pi/3* constants.rho_crit * vals.Omega_m * np.power(8e-3*vals.vk, 3)  # for Vk = 1000km/s, Lp=8Mpc/h
	x = 1 + M_f/M*7.61
	x = np.power(x, -0.526*vals.fd)
	return x