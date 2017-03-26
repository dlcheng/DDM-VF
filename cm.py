import vals
import constants
import numpy as np


def cm_relation(M):
    logm = np.log10(M)
    sigma_logc = 0.132
#	logc_mean = 0.178 - 0.088 * logm  
# Maccio et al. (2008), WMAP5
#	logc_mean = np.log10(cm_planck(M))
    logc_mean = np.log10(cm_bolshoi(M))
    logc = np.random.normal(logc_mean, sigma_logc)
    c = np.power(10, logc)
    # correct for ddm suppression
    c *= cm_ddm_sup(M)
    return c


def cm_planck(M):       # Klypin et al. (2014), Planck
    C0 = 7.4
    gamma = 0.12
    M0 = 5.5e17
    c = C0 * np.power(M / 1e12, -gamma) * (1 + np.power(M / M0, 0.4))
    return c


# M here is the virial mass with Delta = 360
def cm_bolshoi(M):
    return 9.6 * np.power(M / 1e12, -0.075)


def cm_ddm_sup(M):
    M_f = 4.0 * np.pi / 3 * constants.rho_crit * vals.Omega_m * \
        np.power(8e-3 * vals.vk, 3)    # for Vk = 1000km/s, Lp=8Mpc/h
    x = 1 + M_f / M * 20.4
    x = np.power(x, -0.271 * vals.fd)
    return x
