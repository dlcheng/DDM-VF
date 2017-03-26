import vals
import cm
import constants
import numpy as np
from scipy import optimize
from scipy import interpolate


def create_halo_table():
    # how many mass bins
    num = vals.Tb_MF.shape[0]
    vals.Tb_MF[:, 3] = np.array(
        vals.Tb_MF[:, 3] + 0.3, dtype="int64")  # keep the interger part
    table_length = 0
    # maximum sample number
    sample_num = 1000
    for i in range(num):
        table_length += np.min(np.array([vals.Tb_MF[i, 3], sample_num]))
    vals.Tb_halo = np.zeros((table_length, 2))
    pos = 0
    for i in range(num):
        mass = np.power(10, vals.Tb_MF[i, 2])
        loop = np.min(np.array([vals.Tb_MF[i, 3], sample_num]))
        loop = int(loop)
        for j in range(loop):
            # contain random sampling
            c = cm.cm_relation(mass)
            vals.Tb_halo[pos, 0] = max_velocity(mass, c)
            vals.Tb_halo[pos, 1] = vals.Tb_MF[i, 3] / \
                float(loop)       # the weight of this velocity
            pos += 1


def max_velocity(mass, c):
    Delta = 200
    factor = constants.G0 * f_nfw(2.15) / 2.15
    rho = 4.0 * np.pi / 3.0 * Delta * constants.rho_crit * vals.Omega_m
    factor *= np.power(rho, 1.0 / 3.0)
    factor = np.sqrt(factor)
    factor *= 1e-3 * np.sqrt(constants.SUNTOKG / constants.MPCTOM)
    vmax = np.power(mass, 1.0 / 3.0) * np.sqrt(c)
    vmax /= np.sqrt(f_nfw(c))
    # vmax in unit of km/s
    vmax *= factor
    vmax = vals.Sp_infall_cor(vmax)
    return vmax


def f_nfw(x):
    y = np.log(1 + x) - x / (1 + x)
    return y


def create_vf_table(Vmin, Vmax, Num):
    local_bins = np.linspace(np.log10(Vmin), np.log10(Vmax), Num)
    dis = (np.log10(Vmax) - np.log10(Vmin)) / Num
    logv_data = np.log10(vals.Tb_halo[:, 0])
    vals.Hist_VF, vals.Bin_edge = np.histogram(
        logv_data, bins=local_bins, weights=vals.Tb_halo[:, 1])
    vals.Hist_VF /= dis * np.power(vals.boxsize, 3)


# correct Vmax according to Eq. (9) of Klypin et al.(2015)
def baryon_infall_correct(vmax):
    y = vmax / 120.0
    sov = optimize.root(infall_fun, y, (y,), method="hybr")
    return sov.x * 120.0


def infall_fun(x, y):                                             # Eq.(9)
    res = 1 + 0.35 * np.power(x, 6) / (1.0 + np.power(x, 6))
    res = x / res
    return res - y


def init_v_infall_tb():
    vin = np.linspace(1, 500, 500)
    vc = baryon_infall_correct(vin)
    vals.Sp_infall_cor = interpolate.InterpolatedUnivariateSpline(vin, vc)
