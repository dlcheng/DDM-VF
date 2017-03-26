# This routine returns the linear growth factor

import vals
import numpy as np

def g_factor(z):
	x = vals.Omega_m * np.power(1.0+z, 3)
	x += (1.0 - vals.Omega_m - vals.Omega_v) * np.power(1.0+z, 2) + vals.Omega_v
	return x

def unnorm_growth(z):
	g = g_factor(z)
	y1 = vals.Omega_m * np.power(1.0+z, 3)/g
	y2 = vals.Omega_v /g
	x = 1.0 / (1.0 + z) * 2.5 * y1
	x /= np.power(y1, 4.0/7.0) - y2 + (1.0 + 0.5 * y1) * (1.0 + y2/70.0)
	return x

def gr(z):
	return vals.Agr * unnorm_growth(z)