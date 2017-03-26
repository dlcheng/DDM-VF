# All the variables are initialized as float number 0.0. The python mechanism
# allows them to be assigned to the right data structures later in Module
# init.py

# Variables for cosmology
Omega_m = 0.0
Omega_v = 0.0
Omega_b = 0.0
h0 = 0.0

# Varibales for power spectrum
Aps = 0.0                     # the normalization parameter of the PS at z=0
Agr = 0.0                     # the normalization parameter of linear growth
ns = 0.0                      #
sigma8 = 0.0                  #
Ps_flag = 0.0                 # flags for different types of power spectrum
# flag == 1,   scale free
# flag == 2,   BBKS (SCDM, default \Gamma = 0.5)
# flag == 3,   E&H  (power with BAO)
# flag == 4,   Read transfer function from file

# Variables for spline functions used for mass function
Tb_M = 0.0                    # list of the mass of interest
Tb_sigma2 = 0.0               # list of the sigma^2 at z=0 as function of M
# the 1d spline function of log10(sigma2) with log10(M) at z = 0
Sp_logsigma2_logm = 0.0
# the list of transfer function in case of reading from file.
Tb_tk = 0.0
Sp_tk = 0.0                   # the spline funcion of tk from file

# Variables for velocity function for CDM
boxsize = 0.0                 # boxsize of mock halos in unit of Mpc/h
Tb_halo = 0.0                 # the main data structure of halo maximum circular velocity
Tb_MF = 0.0                   # Tb_MF[:,0] is the lower end mass of halo,
# Tb_MF[:,1] is the higher end mass of halo,
# Tb_MF[:,2] is the middle mass of the halo, all mass are in log10 scale
# Tb_MF[:,3] is the number of haloes within the mass bin
Hist_VF = 0.0                 # the histogram of VF
Bin_edge = 0.0                # Bin edage of VF histogram

# Variables for velocity function of DDM
vk = 0.0
fd = 0.0

# Variables for velocity correction of Baryon infall
Sp_infall_cor = 0.0

# Figures
fig = 0.0
ax = 0.0
