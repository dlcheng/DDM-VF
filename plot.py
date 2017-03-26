import vals
import numpy as np
import init
import vf
import matplotlib.pyplot as plt

def vf_klypin(x):
	y = 0.130 * np.power(x/100, -2.9)
	return y

def vf_data(x):
	y = 18.0 / x
	y *= np.exp(-1.0 * np.power(x/250, 3))
	return y

def show_legend(ax, location='best'):
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc=location)	

def plot_bolshoi_vf():
	x = np.linspace(10, 250, 20)
	y = vf_klypin(x)
	vals.ax.plot(x, y, "b-", label="Bolshoi VF")

def plot_obs():
	x = np.linspace(10, 250, 20)
	y = vf_data(x)
	z = y * 1.2
	w = y * 0.8
	vals.fig = plt.figure(figsize=(6, 6))
	vals.ax = vals.fig.add_subplot(1, 1, 1)
	vals.ax.set_xscale("log")
	vals.ax.set_yscale("log")
	vals.ax.set_xlim(10, 210)
	vals.ax.set_ylim(0.03, 10)
	vals.ax.plot(x, y, "k-", label="Klypin et al. (2015)")
	vals.ax.plot(x, z, "k--")
	vals.ax.plot(x, w, "k--")
	show_legend(vals.ax)
    
def plot_theory(vk, fd, Style, Label):
	init.init_ddm_param(vk, fd)
	init.init_halo_mf_table(1e5, 1e16, 100)
	vf.create_halo_table()
	vf.create_vf_table(5, 300, 50)
	vx = np.zeros((vals.Hist_VF.shape[0],))
	dis = vals.Bin_edge[1] - vals.Bin_edge[0]
	for i in range(vx.shape[0]):
		vx[i] = vals.Bin_edge[i] + dis/2
	vx = np.power(10, vx)
	vals.ax.plot(vx, vals.Hist_VF, Style, label=Label)
	show_legend(vals.ax)