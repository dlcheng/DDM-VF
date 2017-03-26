import vals
import init
import growth
import sigma2
import tool
import constants
import tk
import mf
import vf
import cm
import plot
import numpy as np
import matplotlib.pyplot as plt

#reload changed modules

#set parameters (BolshoiP parameters)
#init.set_cosmological_parameters(0.307, 0.693, 0.048, 0.678)
#init.set_power_parameters(0.96, 0.829, 2)

#set parameters (Bolshoi parameters, WMAP7)
init.set_cosmological_parameters(0.27, 0.73, 0.0469, 0.7)
init.set_power_parameters(0.95, 0.82, 4)

#set parameters (Our previous cosmology)
#init.set_cosmological_parameters(0.3, 0.7, 0.047, 0.7)
#init.set_power_parameters(0.96, 0.8, 2)

#set mass function spline
init.init_M_for_MF(1e4, 1e19, 500)

#initialization
init.init_growth_factor()
init.init_power_spectrum()
init.init_logsigma2_slpine()
init.init_boxsize(10000)                # the boxsize can be large to generate enough big haloes,
                                        # while small haloes are studied with subsamples.
vf.init_v_infall_tb()                                        

# Make plots                                         
#plot 1 --- Klypin Eq.(12)
plot.plot_obs()
plot.plot_bolshoi_vf()

#plot 2  --- CDM
plot.plot_theory(100, 0, "r-", "CDM, WMAP7")
#plot.plot_theory(100, 0.5, "r-", r"$V_k=100 km/s, f_d=0.5$") 
#plot.plot_theory(50, 0.5, "g-", r"$V_k=50 km/s, f_d=0.5$")
#plot.plot_theory(20, 0.5, "b-", r"$V_k=20 km/s, f_d=0.5$")
