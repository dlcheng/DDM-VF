# Some tools

import numpy as np
import constants
import vals


def radius(M):
    x = M * 3.0 / 4.0 / np.pi / (constants.rho_crit * vals.Omega_m)
    x = np.power(x, 1.0 / 3.0)
    return x


def mass(R):
    x = 4.0 * np.pi / 3.0 * np.power(R, 3.0) * \
        (constants.rho_crit * vals.Omega_m)
    return x
