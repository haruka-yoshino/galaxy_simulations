###
# THIS IS MAKE_GALAXY WITHOUT UNNECESSARY NOTES AND COMMENTS
# KEEP IT UPDATED FREQUENTLY BY WORKING ONLY IN THE DEV VERSION
# AND COPYING IT OVER
# D O  N O T  WORK HERE 


###########
# IMPORTS #
###########

from my_functions import *
from plotter import plotter
from automate_error_search_etc import *


import numpy as np 

import matplotlib.pyplot as plt 
from matplotlib.ticker import FuncFormatter

from scipy.optimize import fsolve
from scipy.interpolate import interp1d

from tabulate import tabulate

import mpl_scatter_density

import pickle
import os


plt.style.use("/Users/harukayoshino/Documents/uvic/masters1/code/paper_plots.mplstyle")


def main(plot, save_pdf, ID_num, sim_telescope, lookup, bands):

	###################
	# GOOD ISOCHRONES #
	###################

	ALLOWED_m20 = [10.3]

	ALLOWED_m30 = [10.1, 10.15, 10.2, 10.25, 10.3]

	ALLOWED_m35 = [9.75, 9.8, 9.85, 9.9, 9.95, 10.0, 10.05, 10.1, 10.15, 10.2, 10.25, 10.3]

	ALLOWED_m40 = [9.8, 9.85, 9.9, 9.95, 10.0, 10.05, 10.1, 10.15, 10.2, 10.25, 10.3]


	ALLOWED_FEH = [-2.0, -3.0, -3.5, -4.0]

	ALLOWED_AGES = [ALLOWED_m20, ALLOWED_m30, ALLOWED_m35, ALLOWED_m40]


	#############
	# CONSTANTS #
	#############

	D = 16.5*10**6 #np.random.normal(loc = 16.5, scale=2.0)*10**6

	AGE = 10.3

	FEH = -2.0

	# CMD_COLOR_RANGE = [0.8, 1.7]


	####################
	# GALAXY CONSTANTS #
	####################


	filename = "/Users/harukayoshino/Documents/uvic/masters1/ngvs catalogue/Galaxy Catalogue/NGVS_CATALOGUE_FINAL_28Mar2026.csv"

	ID = f"NGVS{ID_num}"

	# NGVS_NAME, RA, DEC, u_mag_tot, g_mag_tot, r_mag_tot, i_mag_tot, z_mag_tot, Re_arcsec, mu_eg, n, PA, ell
	galaxy_parameters = get_parameters(filename, ID)

	# some parameters are used frequently (PA) so they should be a variable like ba and Re


	n = galaxy_parameters["sersic_n"]

	b_n = 2*n - 1/3 + 4/(405*n) + 46 / (25515*n**2) + 131 / (1148175*n**3) - 2194697 / (30690717750*n**4) 
	
	ba = 1 - galaxy_parameters["ell"]

	Re = arcsec_to_pc(galaxy_parameters["Re_arcsec"], D)

	PA_rad = (galaxy_parameters["PA"])*np.pi/180

	data = [["ID", ID, "-"], 
	["Distance", f"{D/(10**6):.1f}", "Mpc"], 
	["Age", f"{10**AGE / (10**9):.1f}", "Gyr"], 
	["[Fe/H]", FEH, "-"], 
	["RA", galaxy_parameters["RA"], "deg"], ["DEC", galaxy_parameters["DEC"], "deg"], 
	["GALFIT", galaxy_parameters["GALFIT"], "-"], ["TYPHON", galaxy_parameters["TYPHON"], "-"], 
	["g", galaxy_parameters["mag_total"]["CFHT_g"], "mag"], 
	["R_e", galaxy_parameters["Re_arcsec"], "arcsec"], ["-", Re, "pc"], 
	["mu_{e,g}", galaxy_parameters["mu_eg"], "mag/arcsec^2"], 
	["n", n, "-"], 
	["PA", galaxy_parameters["PA"], "deg"], ["-", PA_rad, "rad"], 
	["ellipticity", galaxy_parameters["ell"], "-"], 
	["b/a", ba, "-"]]

	headers = ["Parameter", "Value", "Units"]


	print(tabulate(data, headers=headers, tablefmt="grid"))



	#######################
	# TELESCOPE CONSTANTS #
	#######################


	
	ugriz_dictionary = {
	"CFHTugriz": ["CFHT_u", "CFHT_g", "CFHT_r", "CFHT_i_new", "CFHT_z"],
	"SDSS": ["Sloan/U", "Sloan/G", "Sloan/R", "Sloan/I", "Sloan/Z"], 
	"HST_WFC3": ["F336W", "F475W", "F625W", "F775W", "F850LP"], 
	"HST_ACS": [None, "F475W", "F625W", "F775W", "F850LP"]
	}


	IR_dictionary = {
	"HST_WFC3": ["F105W", "F110W", "F125W", "F140W", "F160W"]
	}


	obs_telescope = "HST_WFC3"

	if bands == "ugriz":
		filters = ugriz_dictionary[obs_telescope]
		red = "F775W"
		green = "F475W" 

	elif bands == "IR":
		filters = IR_dictionary[obs_telescope]

		red = "F140W"
		green = "F105W"


	outloc = f"/Users/harukayoshino/Documents/uvic/masters1/simulations/{obs_telescope}/{bands}/"

	
	#################
	# GET ISOCHRONE #
	#################


	isochrone_dictionary = get_isochrone( _age=AGE, _feh=FEH, _telescope=obs_telescope, 
		_filters=filters, _spec_band = bands)#, plot=plot, _xlabel=r"$M_g - M_i$", _ylabel=r"$M_i$") #_IR=False, 

	##################################
	# CONVERT TO APPARENT MAGNITUDES #
	##################################


	m_app = {"D": D}
	# M_abs = {"M_u", "M_g", "M_r", "M_i", "M_z"}
	for label in filters:
		m_app[label] = isochrone_dictionary["M_abs"][label] + 5*np.log10(D) - 5

	################
	# EXTRACT TRGB #
	################


	mag_limit = 28 #np.min(m_app[red])+2 # 28#  
	TRGB = extract_TRGB(m_app[green], m_app[red], mag_limit=mag_limit)#, color_range=CMD_COLOR_RANGE)

	top =  np.argmin(m_app[red][TRGB])
	print(f"TRGB top: {m_app[red][TRGB][top]:.2f}, T_eff = {isochrone_dictionary["T_eff"][TRGB][top]:.2f}")
	bottom = np.argmax(m_app[red][TRGB])
	print(f"TRGB bottom: m_i = {m_app[red][TRGB][bottom]:.2f}, T_eff = {isochrone_dictionary["T_eff"][TRGB][bottom]:.2f}")


	##############
	# GET ERRORS #
	##############


	etc_dictionary = {"HST_WFC3": {
        "ugriz": ["https://etc.stsci.edu/etc/input/wfc3uvis/imaging/", "wfc3_filter_w",
        "sloan"],

        "IR": ["https://etc.stsci.edu/etc/input/wfc3ir/imaging/", "irfilt0",
            "wfc3IR"]
    	}
    }
	min_SNR = 10
	min_t = 90118.744#min_exposure_time(etc_dictionary[obs_telescope][bands], red, min_SNR, isochrone_dictionary["T_eff"][TRGB][bottom] , m_app[red][TRGB][bottom])#_T_eff, _mag)
	print(min_t)
	print(f"minimum exposure time for SNR = 10: {min_t:.2e}s ({min_t/3600:.2e}hrs)")
	

	pklLoc = f'/Users/harukayoshino/Documents/uvic/masters1/pickles/SNR_{min_SNR}_t_{min_t:.2e}'
	if not os.path.exists(pklLoc): os.makedirs(pklLoc)

	pickleFile = pklLoc + f"/D_{D}_maglim_{mag_limit:.2f}_AGE_{AGE}_FEH_{FEH}_{obs_telescope}_{sim_telescope}.pkl"


	# pickleFile = f'/Users/harukayoshino/Documents/uvic/masters1/pickles/SNR_{min_SNR}/D_{D}_maglim_{mag_limit:.2f}_AGE_{AGE}_FEH_{FEH}_{obs_telescope}_{sim_telescope}.pkl'


	if os.path.exists(pickleFile):
		# if it pre-exists, open it
		with open(pickleFile, 'rb') as f:  # Python 3
			min_t, SNR = pickle.load(f)

	else:
		if lookup: 

			exposure_times = [min_t] #[1889901.8655] #[300, 900, 1800, 3600, 5400, 7200, 10800, 18000, 36000] # seconds
			s_list = [2]#[1,2,3,4,5,6,7,8,9,10]
			cols = ["blue", "green", "gold", "orange", "red"]
			exposure_time = exposure_times[-1]
			
			SNR = {}

			for j, t in enumerate(exposure_times):
				print(f"{t:.2e}")
				temp = {}
				for f in filters:#[green, red]: #range(1,5):#enumerate():
					print(f)
					i = ugriz_dictionary[obs_telescope].index(f)
					if ugriz_dictionary[sim_telescope][i] == None:
						pass
					else:
						temp[f] = np.array(error_lookup(etc_dictionary[sim_telescope][bands], 
							ugriz_dictionary[sim_telescope][i], 
							t, isochrone_dictionary["T_eff"][TRGB], m_app[f][TRGB]))#mags_TRGB[i]) #URL, HST_instrument, 

				# print(temp)
				SNR[str(t)] = temp


			with open(pickleFile, 'wb') as f:  # Python 3
				pickle.dump([min_t, SNR], f) # binary file

		else:
			print("Pickle file for errors not found and not looked up")



	################
	# SALPETER IMF #
	################


	dM = np.gradient(isochrone_dictionary["init_mass"])

	SALPETER_NORM = norm_imf(galaxy_parameters["mag_total"]["CFHT_g"], 
		isochrone_dictionary["init_mass"], m_app[green], D, dM, salpeter_imf)

	N_STARS = np.array(SALPETER_NORM*salpeter_imf(isochrone_dictionary["init_mass"])*dM, dtype=int)

	##########################
	# APPLY IMF TO ISOCHRONE #
	##########################


	# OBS_SP_RED, OBS_SP_GREEN, SP_FINAL_MASS 
	sp_dictionary = make_sp(red, green, m_app[red][TRGB], 
		m_app[green][TRGB], 1/SNR[str(min_t)][red], 1/SNR[str(min_t)][green], 
		isochrone_dictionary["final_mass"][TRGB], N_STARS[TRGB])


	
	mag_comp = []
	mag_comp_headers = ["Filter", "NGVS", "Model"]


	for label, val in galaxy_parameters["mag_total"].items():

		i = ugriz_dictionary["CFHTugriz"].index(label)

		obs_label = filters[i]
		mag_comp.append([label, f"{val:.2f}", f"{magnitude(np.sum(N_STARS * flux(m_app[obs_label]))):.2f}"])



	N_TOT = len(sp_dictionary[red])

	mag_comp.append([r"N_{RGB}", "-", f"{N_TOT:.2e}"])
	mag_comp.append([r"N_{tot}", "-", f"{np.sum(N_STARS):.2e}"])


	mask_gt_3Re = (R_ba > 3*Re)

	L_gt_3Re = np.sum(luminosity(sp_dictionary[green][mask_gt_3Re], D))
	L_total = np.sum(luminosity(sp_dictionary[green], D))

	mag_comp.append([r"N_{R>3Re}", "-", f"{len(R_ba[mask_gt_3Re])/len(R_ba):.5f}"])
	mag_comp.append([r"L_{R>3Re}", "-", f"{L_gt_3Re / L_total:.5f}"])


	print(tabulate(mag_comp, headers=mag_comp_headers, tablefmt="grid"))



	R, RMAX = sample_sersic(n, Re, N_TOT, b_n)

	XROT, YROT, R_ba = get_positions(N_TOT, ba, R, PA_rad)

	R_ell = np.array([1,2,3])*Re

	XROT_ell, YROT_ell= get_ellipses(R_ell, ba, PA_rad)

	pos_dictionary = {"R_ba": R_ba, "X": XROT, "Y": YROT, 
			"R_ell": R_ell, "X_ell": XROT_ell, "Y_ell": YROT_ell}



	##################
	# SERSIC PROFILE #
	##################

	# r_input, I_input_fit, r_obs, I_obs, mu_obs
	luminosity_dictionary = get_luminosity_profile(R_ba, ba, n, Re, b_n, RMAX, D, sp_dictionary[green])
	

	if plot:
		square=np.max([np.max(XROT), np.max(YROT)])*1.2 / (10**3)

		plt.figure(figsize=(8,6))
		plt.scatter(XROT/(10**3), YROT/(10**3), s=1, color="black", alpha=0.1)
		# plt.plot(x_major, y_major, 'r-', lw=2)
		# plt.scatter(XROT[~np.isnan(R_ba)][Re_stars], YROT[~np.isnan(R_ba)][Re_stars], s=1, color="red")
		plt.gca().axline((0,0), slope=1/np.tan(PA_rad), color='red', label="semimajor axis")
		
		plt.axis('square')
		plt.title(f"{ID} (n={n})")
		plt.xlabel("x [kpc]")
		plt.ylabel("y [kpc]")
		plt.xlim(-square, square)
		plt.ylim(-square, square)
		plt.gca().invert_xaxis()
		# plt.savefig(outloc)
		plt.show()
		plt.close()


	with open(outloc+f"{ID}_dictionaries.pkl", "wb") as f:
			pickle.dump([D, AGE, FEH, obs_telescope, red, green, mag_limit, pickleFile, 
				galaxy_parameters, sp_dictionary, pos_dictionary, luminosity_dictionary], f)



	if save_pdf:
		plotter(ID, galaxy_parameters["RA"], galaxy_parameters["DEC"], 
			galaxy_parameters["Re_arcsec"], n,
			isochrone_dictionary["M_abs"][green] , isochrone_dictionary["M_abs"][red] , 
			isochrone_dictionary["RGB"] , FEH, AGE, D, 
			SNR[str(min_t)], m_app, filters, isochrone_dictionary["init_mass"], 
			sp_dictionary, TRGB, N_STARS, dM, #sp_dictionary[red], sp_dictionary[green],
			pos_dictionary, PA_rad, Re, luminosity_dictionary, outloc+f"{ID}_sim.pdf")



# for i in range(1,40):
main(plot=False, save_pdf = True, ID_num=1926, sim_telescope="HST_WFC3", lookup = True, bands = "ugriz")


