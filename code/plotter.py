import numpy as np 

import matplotlib.pyplot as plt 
from matplotlib.ticker import FuncFormatter
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import re
import mpl_scatter_density
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d


def plotter(_ID, _RA0, _DEC0, _Re_arcsec, _n, _PA_rad, _Re, _feh, _age, _D, _snr, 
	_filters, _red, _green, 
	_m_app, _isochrone_dictionary, _sp_dictionary, _rgb, _n_stars, _dM, 
	_pos_dictionary, _luminosity_dictionary, _outloc):


	_M_green = _isochrone_dictionary["M_abs"][_green]
	_M_red = _isochrone_dictionary["M_abs"][_red]
	_phase_mask = _isochrone_dictionary["RGB"]
	_init_mass = _isochrone_dictionary["init_mass"]


	_R_ba, _xrot, _yrot, _R_circ, _xrot_circ, _yrot_circ  = _pos_dictionary.values()

	_r_model, _I_model_fit, _r_obs, _I_obs, _mu_obs = _luminosity_dictionary.values()

	_obs_sp_red, _obs_sp_green, _sp_final_mass = _sp_dictionary.values()

	 
	# _sp_final_mass, _rgb, _n_stars, _dM, 
	# _obs_sp_red, _obs_sp_green, _R_ba, _xrot, _yrot, _R_circ, _xrot_circ, _yrot_circ,  #_x_circ, _y_circ,
	# _PA_rad, _Re, _r_model, _I_model_fit, _r_obs, _I_obs, _mu_obs, 




	#Set up the plots
	gs=gridspec.GridSpec(3, 8, height_ratios=[1, 0.1, 1], width_ratios=[1, 0.4, 1, 0.4, 1, 0.6, 1, 0.2])
	gs.update(left=0.1, right=0.98, bottom=0.1, top=0.90, wspace=0.2, hspace=0.3)

	

	fig = plt.figure(figsize=(20,9))

	ax11 = fig.add_subplot(gs[0,0]) 
	ax12 = fig.add_subplot(gs[0,2]) 
	ax13 = fig.add_subplot(gs[0,4]) 
	ax14 = fig.add_subplot(gs[0,6])
	ax21 = fig.add_subplot(gs[2,0]) 
	ax22 = fig.add_subplot(gs[2,2], projection='scatter_density')
	ax23 = fig.add_subplot(gs[2,4], projection='scatter_density')
	ax24 = fig.add_subplot(gs[2,6])

	fig.suptitle(_ID, fontsize=20)





	ax11.scatter(_M_green - _M_red, _M_red, s=5, color="blue", label="isochrone")
	ax11.scatter(_M_green[_phase_mask] - _M_red[_phase_mask], _M_red[_phase_mask], s=5, color="red", label="RGB")
	ax11.set_xlabel(r"$M_g - M_i$")
	ax11.set_ylabel(r"$M_i$")
	# ax11.legend()
	ax11.set_title(f"[Fe/H] = {_feh}\nlog_10(Age) = {_age}")
	ax11.invert_yaxis()




	colors = np.where(_phase_mask[_rgb], "red", "blue")
	ax12.scatter(np.log10(_init_mass[_rgb]), np.log10(_n_stars[_rgb]/_dM[_rgb]), s=5, c=colors)


	xticks = np.linspace(np.min(_init_mass[_rgb]), np.max(_init_mass[_rgb]), 3)
	ax12.set_xticks(np.log10(xticks), labels=[f'{x:.4f}' for x in xticks])
	ax12.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{10**x:.6f}"))




	ax12.set_xlabel(r"$M_{*, initial}~[M_{\odot}]$")
	ax12.set_ylabel(r"$\log_{10}(N/\Delta M)$")
	ax12.set_title("Salpeter IMF (red => RGB)")




	ax13_plot = ax13.scatter(_obs_sp_green - _obs_sp_red, _obs_sp_red, s=5, c=_sp_final_mass, cmap="viridis")

	ax13c = inset_axes(ax13, width='5%', height='100%', loc='lower left', bbox_to_anchor=(1.05, 0, 1, 1), bbox_transform=ax13.transAxes, borderpad=0)
	ax13cb = fig.colorbar(ax13_plot, cax=ax13c, label=r"$M_{*, current}~[M_{\odot}]$")

	ax13.set_xlabel(r"$m_g - m_i$")
	ax13.set_ylabel(r"$m_i$")
	ax13.set_title(r"Simulated CMD")
	ax13.invert_yaxis()






	cols = ["blue", "green", "gold", "orange", "red"]
	for _i, _f in enumerate(_filters):
		ax14.scatter(_m_app[_f][_rgb], _snr[_f], s=2, color = cols[_i], label=_f)


	
	# for t, dic in _snr.items():
	# 	j=0
	# 	for f, val in dic.items():
	# 		ax14.scatter(_m_app[f][_rgb], dic[f], s=2, color=cols[j])
	# 		j+=1


	ax14.set_title("SNR vs magnitude")
	ax14.set_xlabel("apparent magnitude")
	ax14.set_ylabel("SNR")
	ax14.legend()



	# to_primary = interp1d(_mu_obs, _I_obs, fill_value="extrapolate")
	# to_secondary = interp1d(_I_obs, _mu_obs,  fill_value="extrapolate")


	ax21.plot(_r_model, _I_model_fit, color="black", label="theoretical profile")
	ax21.errorbar(_r_obs, _I_obs, yerr=np.sqrt(10**_I_obs), fmt='o', linestyle='none', color="red", label="binned profile")
	ax21.vlines(x=np.log10(_Re), ymin=np.min(_I_obs), ymax = np.max(_I_obs), label=r"$mR_e$")
	ax21.set_title(f"Sersic profile (n = {_n})")
	ax21.set_xlabel(r"$\log_{10}(R)$ [pc]")
	ax21.set_ylabel(r"$\log_{10}$(number / area) [#/pc$^2$]")  



	ax21sec = ax21.twinx()
	ax21sec.set_ylim(ax21.get_ylim())

	# choose fewer evenly spaced ticks
	yticks = np.linspace(ax21.get_ylim()[0], ax21.get_ylim()[1], 5)

	idx = np.argsort(_I_obs)
	labels = np.interp(yticks, _I_obs[idx], _mu_obs[idx])

	ax21sec.set_yticks(yticks)
	ax21sec.set_yticklabels([f"{v:.2f}" for v in labels])


	# secax_y = ax21.secondary_yaxis('right', functions=(to_secondary, to_primary))
	ax21sec.set_ylabel(r"$\mu_{e}$ [mag / arcsec$^2$]")  




	

	
	_square = np.max([-np.min(_xrot), np.max(_xrot), -np.min(_yrot), np.max(_yrot)]) / _Re

	# bin_width = 10**2 / _Re # 100 pc resolution 
	# bins_x = np.arange(-_square, _square + bin_width, bin_width)
	# bins_y = np.arange(-_square, _square + bin_width, bin_width)

	# print(len(bins_x), len(bins_y))


	# H, xedges, yedges = np.histogram2d(_xrot/_Re, _yrot/_Re, bins = [bins_x, bins_y])

	# dx = xedges[1] - xedges[0]
	# dy = yedges[1] - yedges[0]
	# density = H.T / (dx * dy)

	# xcenters = (xedges[:-1] + xedges[1:]) / 2
	# ycenters = (yedges[:-1] + yedges[1:]) / 2


	# density23 = ax22.imshow(density, origin='lower', 
	# 	extent = [-_square, _square, -_square, _square], #interpolation='nearest',  

	# 	aspect="equal", cmap = "viridis", norm="log")
	
	density22 = ax22.scatter_density(_xrot/_Re, _yrot/_Re, cmap="viridis")

	fig.canvas.draw()

	vals22 = density22.get_array()
	vals22 = vals22[np.isfinite(vals22) & (vals22 > 0)]

	vmin22 = np.percentile(vals22, 0.01)
	vmax22 = np.percentile(vals22, 99.99)
	density22.set_norm(LogNorm(vmin=vmin22, vmax=vmax22))



	ax22.scatter(_xrot_circ/_Re, _yrot_circ/_Re, s=1, label=r"ellipses at $mR_e$", color="red")
	ax22.axline((0,0), slope=1/np.tan(_PA_rad), color='black', label="semimajor axis")
	ax22.set_xlim(-_square, _square)
	ax22.set_ylim(-_square, _square)
	ax22.set_xlabel(r"$x ~[R_e]$")
	ax22.set_ylabel(r"$y ~[R_e]$")
	ax22.set_title("2D distribution of stars")
	# ax[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))



	secax_x = ax22.secondary_xaxis('top',
	                             functions=(lambda x: _Re_arcsec / 3600*x + _RA0, lambda x: (x - _RA0) / _Re_arcsec / 3600))#pc_to_arcsec, arcsec_to_pc))

	secax_y = ax22.secondary_yaxis('right',
	                             functions=(lambda x: _Re_arcsec / 3600*x + _DEC0, lambda x: (x - _DEC0) / _Re_arcsec / 3600))

	secax_x.set_xlabel("RA (deg)")
	secax_y.set_ylabel("DEC (deg)")

	ax22.invert_xaxis()



	
	_square = np.max([-np.min(_xrot), np.max(_xrot), -np.min(_yrot), np.max(_yrot)]) / (10**3)


	cmap = plt.get_cmap('viridis')
	darkest_color = cmap(0.0)
	ax23.set_facecolor(darkest_color)

	# bin_width = 10**2 / (10**3) # 100 pc resolution, 1.25 arcsec resolution
	# bins_x = np.arange(-_square, _square + bin_width, bin_width)
	# bins_y = np.arange(-_square, _square + bin_width, bin_width)


	# H, xedges, yedges = np.histogram2d(_xrot/(10**3), _yrot/(10**3), bins = [bins_x, bins_y])

	# dx = xedges[1] - xedges[0]
	# dy = yedges[1] - yedges[0]
	# density = H.T / (dx * dy)

	# xcenters = (xedges[:-1] + xedges[1:]) / 2
	# ycenters = (yedges[:-1] + yedges[1:]) / 2


	# density23 = ax23.imshow(density, origin='lower', 
	# 	extent = [-_square, _square, -_square, _square],  #interpolation='nearest',  
	# 	cmap = "viridis", norm="log", aspect="equal")

	

	density23 = ax23.scatter_density(_xrot / (10**3), _yrot / (10**3), cmap="viridis")

	fig.canvas.draw()

	vals23 = density23.get_array()
	vals23 = vals23[np.isfinite(vals23) & (vals23 > 0)]

	vmin = np.percentile(vals23, 0.01)
	vmax = np.percentile(vals23, 99.99)
	density23.set_norm(LogNorm(vmin=vmin, vmax=vmax))

	ax23c = inset_axes(ax23, width='5%', height='100%', loc='lower left', bbox_to_anchor=(1.4, 0, 1, 1), bbox_transform=ax23.transAxes, borderpad=0)
	ax23cb = fig.colorbar(density23, cax=ax23c, label="number density")

	ax23.set_xlim(-_square, _square)
	ax23.set_ylim(-_square, _square)
	ax23.set_xlabel(r"$x ~[kpc]$")
	ax23.set_ylabel(r"$y ~[kpc]$")
	ax23.set_title("2D distribution of stars")



	secax_x = ax23.secondary_xaxis('top',
	                             functions=(lambda x: np.arctan(x*(10**3) / _D) * 3600 * 180 / np.pi, lambda x: _D / (10**3) * np.tan(x* np.pi / 3600 / 180)))#pc_to_arcsec, arcsec_to_pc))

	secax_y = ax23.secondary_yaxis('right',
	                             functions=(lambda x: np.arctan(x*(10**3) / _D) * 3600 * 180 / np.pi, lambda x: _D / (10**3) * np.tan(x* np.pi / 3600 / 180)))

	secax_x.set_xlabel("arcseconds")
	secax_y.set_ylabel("arcseconds")
	ax23.invert_xaxis()


	ax24.plot([0,1], [0,1], color="black")
	ax24.plot([1,0], [0,1], color="black")
	ax24.set_xticks([])
	ax24.set_yticks([])
	plt.savefig(_outloc)

	plt.close()

