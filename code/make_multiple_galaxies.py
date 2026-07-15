
# THIS IS THE DEV VERSION
# DO ALL THE WORK HERE
# COPY IT OVER TO THE COMP VERSION WHEN IT WORKS


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
from scipy.special import gammaincinv

from tabulate import tabulate

import mpl_scatter_density

import pickle
import os
import time



plt.style.use("/Users/harukayoshino/Documents/uvic/masters1/code/paper_plots.mplstyle")


##########################
# START OF MAIN FUNCTION #
##########################


def main(plot, save_pdf, lookup, bands, magrange, user_min_t):

	#############
	# CONSTANTS #
	#############

	D = 16.5*10**6 

	AGE = 10.3

	FEH = -2.0


	#######################
	# TELESCOPE CONSTANTS #
	#######################

	
	ugriz_dictionary = {
	"CFHTugriz": ["CFHT_u", "CFHT_g", "CFHT_r", "CFHT_i_new", "CFHT_z"],
	"SDSS": ["Sloan/U", "Sloan/G", "Sloan/R", "Sloan/I", "Sloan/Z"], 
	"HST_WFC3": ["F336W", "F475W", "F625W", "F775W", "F850LP"]#, 
	# "HST_ACS": [None, "F475W", "F625W", "F775W", "F850LP"]
	}


	IR_dictionary = {
	"HST_WFC3": ["F105W", "F110W", "F125W", "F140W", "F160W"], 
	"WFIRST": ["R062", "Z087", "Y106", "J129",  "H158" ]#, "F184", ] #"W146",
	}


	# obs_telescope = "HST_WFC3"
	obs_telescope = "WFIRST"

	if bands == "ugriz":
		filters = ugriz_dictionary[obs_telescope]
		red = "F775W" 
		green = "F475W" 

	elif bands == "IR":
		filters = IR_dictionary[obs_telescope]

		# red = "F160W" 
		# green = "F125W" 


		red = "H158" #"F140W" ##"CFHT_i_new" 
		green = "J129" #"F105W" #"CFHT_g" 

	mag_limit = 27.52 #28 #np.min(m_app[red])+2 #28#  

	#################
	# GET ISOCHRONE #
	#################

	# this might need to be done multiple times for all the isochrones that we need to make the galaxies
	# and then we can choose the one we need in the galaxy loop
	isochrone_dictionary = get_isochrone( _age=AGE, _feh=FEH, _telescope=obs_telescope, 
		_filters=filters, _spec_band = bands)#, plot=plot, _xlabel=r"$M_g - M_i$", _ylabel=r"$M_i$") #_IR=False, 

	CFHT_M_abs = get_CFHT_mags( _age=AGE, _feh=FEH, _telescope="CFHTugriz", 
		_filters=ugriz_dictionary["CFHTugriz"], _spec_band = "ugriz")


	##################################
	# CONVERT TO APPARENT MAGNITUDES #
	##################################

	m_app = {}
	# M_abs = {"M_u", "M_g", "M_r", "M_i", "M_z"}
	for label in isochrone_dictionary["M_abs"].keys():
		m_app[label] = isochrone_dictionary["M_abs"][label] + 5*np.log10(D) - 5


	CFHT_m_app = {}
	for label in CFHT_M_abs.keys():
		CFHT_m_app[label] = CFHT_M_abs[label] + 5*np.log10(D) - 5

	################
	# EXTRACT TRGB #
	################


	
	TRGB = extract_TRGB(m_app[green], m_app[red], mag_limit=mag_limit)#, color_range=CMD_COLOR_RANGE)

	top =  np.argmin(m_app[red][TRGB])
	print(f"TRGB top: {m_app[red][TRGB][top]:.2f}, T_eff = {isochrone_dictionary["T_eff"][TRGB][top]:.2f}")
	bottom = np.argmax(m_app[red][TRGB])
	print(f"TRGB bottom: m_i = {m_app[red][TRGB][bottom]:.2f}, T_eff = {isochrone_dictionary["T_eff"][TRGB][bottom]:.2f}")


	##############
	# GET ERRORS #
	##############

	
	# 90118.744 #
	min_SNR = "_"


	# min_t = 90118.744 #min_exposure_time(etc_dictionary[obs_telescope][bands], red, min_SNR, isochrone_dictionary["T_eff"][TRGB][bottom] , m_app[red][TRGB][bottom])#_T_eff, _mag)
	# print(f"minimum exposure time for SNR = {min_SNR}: {min_t:.2e}s ({min_t/3600:.2e}hrs)")
	# exposure_times = [min_t] 
	# s_list = [2]
	# cols = ["blue", "green", "gold", "orange", "red"]
	# exposure_time = exposure_times[-1]
	# SNR = {}

	if obs_telescope == "HST_WFC3":
		etc_dictionary = {"HST_WFC3": {
	        "ugriz": ["https://etc.stsci.edu/etc/input/wfc3uvis/imaging/", "wfc3_filter_w",
	        "sloan"],

	        "IR": ["https://etc.stsci.edu/etc/input/wfc3ir/imaging/", "irfilt0",
	            "wfc3IR"]
	    	}
	    }

		# this needs some improvement
		pklLoc = f'/Users/harukayoshino/Documents/uvic/masters1/pickles/SNR_{min_SNR}'
		if not os.path.exists(pklLoc): os.makedirs(pklLoc)

		# this should include the isochrone that is being simulated
		pickleFile = pklLoc + f"/pls_use_this_{bands}_1hr.pkl"


		if os.path.exists(pickleFile):
			# if it pre-exists, open it
			with open(pickleFile, 'rb') as f:  # Python 3
				min_t, SNR = pickle.load(f)


		else:
			if lookup: 
				min_t = {}

				SNR = {}

				# for each filter,
				for i, f in enumerate(filters):

					# if the minimum exposure time is not specified already, 
					if user_min_t[i] == 0:

						# get the minimum exposure time needed to achieve a minimum SNR
						min_t[f] = min_exposure_time(etc_dictionary[obs_telescope][bands], f, min_SNR, 
							isochrone_dictionary["T_eff"][TRGB][bottom] , m_app[f][TRGB][bottom])

					else:
						min_t[f] = user_min_t[i]


					SNR[f] = np.array(error_lookup(etc_dictionary[obs_telescope][bands], 
						f, min_t[f], isochrone_dictionary["T_eff"][TRGB], m_app[f][TRGB]))

				with open(pickleFile, 'wb') as f:  # Python 3
					pickle.dump([min_t, SNR], f) # binary file

			else:
				print("Pickle file for errors not found and not looked up")

				exit()
		# for f in filters:
		# 	plt.scatter(m_app[f][TRGB], SNR[f], label=f)
		# 	print(np.max(SNR[f]))

		# plt.xlabel("apparent magnitude")
		# plt.ylabel("SNR")
		# plt.title("minimum SNR = 10 in all bands")
		# plt.legend()
		# # plt.show()
		# plt.close()

		# # exit()

		# plt.figure(figsize=(9,6))
		# for f in filters:
		# 	plt.scatter(f, min_t[f] / 3600)
		# 	print(min_t[f]/3600)
		# plt.xlabel(f"{bands} filter")
		# plt.ylabel("exposure time (hours)")
		# plt.title("minimum exposure time needed for SNR = 10")
		# # plt.show()
		# plt.close()


	
	else:
		pickleFile = "linear"
		SNR = {}

		# 1 hour
		# interp_red = interp1d([26.86, 27.47], [7.42, 4.27], kind='linear', fill_value="extrapolate")

		# 10 hours
		interp_red = interp1d([26.86, 27.47], [23.69, 13.63], kind='linear', fill_value="extrapolate")
		SNR[red]= interp_red(m_app[red][TRGB]) 

		# 1 hour
		# interp_green = interp1d([27.04, 27.62], [5.87, 3.43], kind='linear', fill_value="extrapolate")

		# 10 hours
		interp_green = interp1d([27.04, 27.62], [18.75, 10.95], kind='linear', fill_value="extrapolate")
		SNR[green] = interp_green(m_app[green][TRGB]) # what it should actually be


		# pickleFile = "uniform"
		# SNR[f] = np.ones(len(m_app[red][TRGB]))*5

	####################
	# GALAXY CONSTANTS #
	####################

	
	# exit()

	filename = "/Users/harukayoshino/Documents/uvic/masters1/ngvs catalogue/Galaxy Catalogue/NGVS_CATALOGUE_FINAL_28Mar2026.csv"



	names = get_detections(filename, magrange)

	print(len(names))
	for ID in names:

		print(f"\n{ID}\n")

		# NGVS_NAME, RA, DEC, u_mag_tot, g_mag_tot, r_mag_tot, i_mag_tot, z_mag_tot, Re_arcsec, mu_eg, n, PA, ell

		outloc = f"/Users/harukayoshino/Documents/uvic/masters1/simulations/{obs_telescope}/maglim_{mag_limit}/{bands}/{green}_{red}/exp_time_10hr/{ID}"

		if not os.path.exists(outloc): os.makedirs(outloc)


		### 
		# since we already extracted necessary shit before the loop we could just replace this with indices n shit
		###

		galaxy_parameters = get_parameters(filename, ID)

		# some parameters are used frequently (PA) so they should be a variable like ba and Re

		n = galaxy_parameters["sersic_n"]
		
		b_n = 2*n - 1/3 + 4/(405*n) + 46 / (25515*n**2) + 131 / (1148175*n**3) - 2194697 / (30690717750*n**4)

		ba = 1 - galaxy_parameters["ell"]

		Re = arcsec_to_pc(galaxy_parameters["Re_arcsec"], D)

		PA_rad = (galaxy_parameters["PA"])*np.pi/180

		R_MAX = Re*(np.log(10**-10) / (-b_n) + 1) ** n
		if np.isnan(R_MAX):
			try_b_n = gammaincinv(2*n, 0.5)
			try_R_MAX = Re*(np.log(10**-10) / (-try_b_n) + 1) ** n

			if np.isnan(try_R_MAX):
				print("b_n is incompatible")
				continue
			else:
				b_n = try_b_n
				R_MAX = try_R_MAX

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

		with open(outloc+f"/{ID}_params.txt", "w") as f:
		    f.write(tabulate(data, headers=headers, tablefmt="grid",  disable_numparse=True))

		
		# LOOP THROUGH LIST 

		################
		# SALPETER IMF #
		################


		dM = np.gradient(isochrone_dictionary["init_mass"])

		SALPETER_NORM = norm_imf(galaxy_parameters["mag_total"]["CFHT_g"], 
			isochrone_dictionary["init_mass"], CFHT_m_app["CFHT_g"], D, dM, salpeter_imf)

		N_STARS = np.array(SALPETER_NORM*salpeter_imf(isochrone_dictionary["init_mass"])*dM, dtype=int)


		if plot:
			print("insert updated plots here")


		##########################
		# APPLY IMF TO ISOCHRONE #
		##########################

		sp_dictionary = make_sp(red, green, m_app[red][TRGB], 
			m_app[green][TRGB], 1/SNR[red], 1/SNR[green], 
			isochrone_dictionary["final_mass"][TRGB], N_STARS[TRGB])

		
		N_TOT = np.sum(N_STARS[TRGB]) #len(sp_dictionary[red])


		
		R = sample_sersic(n, Re, N_TOT, b_n, R_MAX)

		XROT, YROT, R_ba = get_positions(N_TOT, ba, R, PA_rad)

		R_ell = np.array([1,2,3])*Re

		XROT_ell, YROT_ell= get_ellipses(R_ell, ba, PA_rad)

		pos_dictionary = {"R_ba": R_ba, "X": XROT, "Y": YROT, 
				"R_ell": R_ell, "X_ell": XROT_ell, "Y_ell": YROT_ell}



		mag_comp = []
		mag_comp_headers = ["Filter", "NGVS", "Model"]


		for label, val in galaxy_parameters["mag_total"].items():

			mag_comp.append([label, f"{val:.2f}", f"{magnitude(np.sum(N_STARS * flux(CFHT_m_app[label]))):.2f}"])


		for f in filters:
			mag_comp.append([f, "-", f"{magnitude(np.sum(N_STARS * flux(m_app[f]))):.2f}"])




		mag_comp.append([r"N_{RGB}", "-", f"{N_TOT:.2e}"])
		mag_comp.append([r"N_{tot}", "-", f"{np.sum(N_STARS):.2e}"])


		mask_gt_3Re = (R_ba > 3*Re)

		L_gt_3Re = np.sum(luminosity(sp_dictionary[green][mask_gt_3Re], D))
		L_total = np.sum(luminosity(sp_dictionary[green], D))

		mag_comp.append([r"N_{R>3Re}", "-", f"{len(R_ba[mask_gt_3Re])/len(R_ba):.5f}"])
		mag_comp.append([r"L_{R>3Re}", "-", f"{L_gt_3Re / L_total:.5f}"])

		with open(outloc+f"/{ID}_model_mag_tot.txt", "w") as f:
		    f.write(tabulate(mag_comp, headers=mag_comp_headers, tablefmt="grid", disable_numparse=True))


		if plot:
			print("insert updated plots here")


		##################
		# SERSIC PROFILE #
		##################

		
		
		luminosity_dictionary = get_luminosity_profile(R_ba, ba, n, Re, b_n, R_MAX, D, sp_dictionary[green])
			

		# square = np.max([np.max(XROT), np.max(YROT)]) +5000
		square = np.max([-np.min(XROT), np.max(XROT), -np.min(YROT), np.max(YROT)]) / (10**3)

		plt.figure(figsize=(8,6))
		plt.scatter(XROT/(10**3), YROT/(10**3), s=1, color="black", alpha=0.1)
		plt.gca().axline((0,0), slope=1/np.tan(PA_rad), color='red', label="semimajor axis")
		
		plt.axis('square')
		plt.title(f"{ID} (n = {n}, PA = {galaxy_parameters["PA"]})")
		plt.xlabel("x [kpc]")
		plt.ylabel("y [kpc]")
		plt.xlim(-square, square)
		plt.ylim(-square, square)

		plt.gca().invert_xaxis()
		plt.savefig(outloc+f"/{ID}_2D.png")
		# plt.show()
		plt.close()

		# exit()

		if plot:
			print("insert updated plots here")


		with open(outloc+f"/{ID}_locs.pkl", 'wb') as f:  # Python 3
			pickle.dump([ID, galaxy_parameters["RA"], galaxy_parameters["DEC"], 
				galaxy_parameters["Re_arcsec"], 
				XROT, YROT, Re, ba, PA_rad], f) # binary file

		with open(outloc+f"/{ID}_dictionaries.pkl", "wb") as f:
			pickle.dump([D, AGE, FEH, obs_telescope, red, green, mag_limit, pickleFile, 
				galaxy_parameters, sp_dictionary, pos_dictionary, luminosity_dictionary], f)

		# save all of these parameters in a good way somehow
		# some of them are probably unnecessary
		# others can be pushed through in dictionaries and opened later
		# diagnostic_plots should be {ID}_diagnostic_plots

		if save_pdf or ID == "NGVS1362":
			plotter(ID, galaxy_parameters["RA"], galaxy_parameters["DEC"], 
				galaxy_parameters["Re_arcsec"], n, PA_rad, Re, FEH, AGE, D, SNR,
				filters, red, green, m_app, isochrone_dictionary, sp_dictionary, 
				TRGB, N_STARS, dM, pos_dictionary, luminosity_dictionary, 
				outloc+f"/{ID}_diagnostic_plots.png")
			
			# plotter(ID, galaxy_parameters["RA"], galaxy_parameters["DEC"], 
			# 	galaxy_parameters["Re_arcsec"], n, PA_rad, Re, FEH, AGE, D, SNR,
			# 	[red, green], red, green, m_app, isochrone_dictionary, sp_dictionary, 
			# 	TRGB, N_STARS, dM, pos_dictionary, luminosity_dictionary, 
			# 	outloc+f"/{ID}_diagnostic_plots.png")

# for i in range(1,40):
start_time = time.time()

main(plot=False, save_pdf = False, lookup = True, bands = "IR", magrange = [8,22], user_min_t = [3600, 3600, 3600, 3600, 3600]) #ID_num=1362,

print(f"--- {(time.time() - start_time)} seconds ---")



