# This routine calculates the sigma^2 as function of M and a directly from
# its definition

import vals
import tk
import tool
import numpy as np
from scipy import integrate

# Sigma^2 calculated at z=0


def calculate_sigma2(M):
    R = tool.radius(M)
    return vals.Aps * unnorm_sigma2(R)

# The Fourier transform of top-hat window


def win_th(y):
    x = 1
    if y <= 1e-8:
        x = 1.0 - 0.1 * y * y
    else:
        x = 3.0 / np.power(y, 3) * (np.sin(y) - y * np.cos(y))
    return x


def sigma2_int_kernel(k, R):
    x = 0.5 / np.pi / np.pi
    x *= np.power(k, vals.ns + 2.0)
    x *= np.power(tk.prim_tk(k), 2)
    x *= np.power(win_th(k * R), 2)
    return x


def unnorm_sigma2(R):
    x, err = integrate.quad(sigma2_int_kernel, 0, np.inf,
                            (R,), epsabs=0.0, epsrel=1e-5)
    return x


def norm_ps():
    if vals.Ps_flag == 1:
        vals.Aps = 1.0
    else:
        vals.Aps = vals.sigma8 * vals.sigma8 / unnorm_sigma2(8.0)
