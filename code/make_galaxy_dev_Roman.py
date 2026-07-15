

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


plt.style.use("/Users/harukayoshino/Documents/uvic/masters1/code/paper_plots.mplstyle")

#############
# FUNCTIONS #
#############


##########################
# START OF MAIN FUNCTION #
##########################

"""
Ok so right now we have two options for the AGE and FEH variables
1. they are pre-defined (they need to be changed in the code before running)
2. they are defined by the user (the code doesn't need to be changed before running)

The problem is that we don't know which combination of AGE and FEH will give a good 
isochrone (clear TRGB).

But this is actually completely fine, because later, we need to be able to make a 
model with a bunch of combinations of ages and metallicites and compare it to 
observational data. This means that the age and metallicity would be parameters in 
a MCMC (if it were to be done this way).

So for now, what we can do is manually find the combination of ages and metallicities
that give us a good TRGB CMD, and only allow those. Or, even simpler, we can just pre-
define one that we like.
"""


def main(plot, save_pdf, ID_num, lookup, bands, user_min_t):

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




###
# MAKE A DICTIONARY FOR THE GALAXY PARAMETERS
# DO NOT INCLUDE DISTANCE, AGE, AND FE/H (THOSE ARE IN STELLAR_POPULATION AND ISOCHRONE) 
###
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

	R_MAX = Re*(np.log(10**-10) / (-b_n) + 1) ** n
	if np.isnan(R_MAX):
		try_b_n = gammaincinv(2*n, 0.5)
		try_R_MAX = Re*(np.log(10**-10) / (-try_b_n) + 1) ** n

		if np.isnan(try_R_MAX):
			print("b_n is incompatible")
			exit()
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


	print(tabulate(data, headers=headers, tablefmt="grid"))



	#######################
	# TELESCOPE CONSTANTS #
	#######################

####
# MAKE A DICTIONARY FOR THE FILTERS
# THE ETC-SPECIFIC NAMES CAN BE IGNORED FOR NOW, IT'S A SEPARATE VARIABLE
####


	
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



	# telescope that did the observations
	# obs_telescope="CFHTugriz"
	# obs_telescope = "HST_WFC3"
	obs_telescope = "WFIRST"
	# as opposed to the telescope / instrument that's being simulated
	# should make that distinction clearer

	if bands == "ugriz":
		filters = ugriz_dictionary[obs_telescope]
		red = "F775W"  ##"CFHT_i_new" #
		green = "F475W"   #"CFHT_g" #

	elif bands == "IR":
		filters = IR_dictionary[obs_telescope]

		red = "H158" #"F140W" ##"CFHT_i_new" 
		green = "J129" #"F105W" #"CFHT_g" 

	# filters = ["CFHT_u", "CFHT_g", "CFHT_r", "CFHT_i_new", "CFHT_z"]
	# red_band="CFHT_i_new"
	# green_band="CFHT_g"



	
	#################
	# GET ISOCHRONE #
	#################

	isochrone_dictionary = get_isochrone( _age=AGE, _feh=FEH, _telescope=obs_telescope, 
		_filters=filters, _spec_band = bands)#, plot=plot, _xlabel=r"$M_g - M_i$", _ylabel=r"$M_i$") #_IR=False, 

	CFHT_M_abs = get_CFHT_mags( _age=AGE, _feh=FEH, _telescope="CFHTugriz", 
		_filters=ugriz_dictionary["CFHTugriz"], _spec_band = "ugriz")



	# some variables are used frequently (initial mass, effective temperature) so they should be variables

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
	# SALPETER IMF #
	################


	dM = np.gradient(isochrone_dictionary["init_mass"])

	SALPETER_NORM = norm_imf(galaxy_parameters["mag_total"]["CFHT_g"], 
		isochrone_dictionary["init_mass"], CFHT_m_app["CFHT_g"], D, dM, salpeter_imf)

	N_STARS = np.array(SALPETER_NORM*salpeter_imf(isochrone_dictionary["init_mass"])*dM, dtype=int)

	################
	# EXTRACT TRGB #
	################
	"""
	Do we do this by color and magnitude cuts, or do we use the "phase" mask that's given in the
	isochrone?
	"""

	mag_limit = 27.52 #np.min(m_app[red])+2 # 28#  

	TRGB = extract_TRGB(m_app[green], m_app[red], mag_limit=mag_limit)#, color_range=CMD_COLOR_RANGE)

	top =  np.argmin(m_app[red][TRGB])
	print(f"TRGB top: {m_app[red][TRGB][top]:.2f}, T_eff = {isochrone_dictionary["T_eff"][TRGB][top]:.2f}")
	print(f"TRGB top: {m_app[green][TRGB][top]:.2f}")#, T_eff = {isochrone_dictionary["T_eff"][TRGB][top]:.2f}")
	bottom = np.argmax(m_app[red][TRGB])
	print(f"TRGB bottom: m_i = {m_app[red][TRGB][bottom]:.2f}, T_eff = {isochrone_dictionary["T_eff"][TRGB][bottom]:.2f}")
	print(f"TRGB top: {m_app[green][TRGB][bottom]:.2f}")#
	exit()

###
# MAKE A DICTIONARY FOR THE RGB STARS (AFTER EXTRACTING THE RGB STARS, BUT BEFORE APPLYING THE IMF)
# SHOULD BE JUST 10 OR SO UNIQUE STARS AT THE TOP OF THE ISOCHRONE
# SOME KEYS: 'INIT_MASS', 'FINAL_MASS', 'T_EFF', 'M_APP_U', 'SNR_U', 'ERR_U', 'N_STARS'
###
	##############
	# GET ERRORS #
	##############


	###
	# LOOK UP ERROR IN ERROR TABLE
	# BUT ALSO RUN AUTOMATE ERROR SEARCH AND GET THE ERROR FOR EACH TEMPERATURE AND MAGNITUDE
	###


	# etc_dictionary = {"HST_WFC3": ["https://etc.stsci.edu/etc/input/wfc3uvis/imaging/", "wfc3_filter_w"],
	# "HST_ACS": ["https://etc.stsci.edu/etc/input/acs/imaging/", "wfcfilt0"]} 


	etc_dictionary = {"HST_WFC3": {
        "ugriz": ["https://etc.stsci.edu/etc/input/wfc3uvis/imaging/", "wfc3_filter_w",
        "sloan"],

        "IR": ["https://etc.stsci.edu/etc/input/wfc3ir/imaging/", "irfilt0",
            "wfc3IR"]
    	}
    }

	###
	# IF WE USE HST WFC3 ETC WITH CFHT ISOCHRONE (AGE 10.3 FEH -2.0 MAGLIM = TRGB+2),
	# THE TIME NEEDED FOR THE MAGLIM TO REACH SNR = 10 IS 1,889,901.8655 SECONDS (21 DAYS)

	# SAME THING, BUT IF WE USE THE CORRESPONDING HST_WFC ISOCHRONE WITH F775W, 
	# TIME IS 886,332.8758 SECONDS (10 DAYS)
	###
	# 90118.744 #


	min_SNR = 10
	# min_t = min_exposure_time(etc_dictionary[obs_telescope][bands], red, min_SNR, 
	# 			isochrone_dictionary["T_eff"][TRGB][bottom] , m_app[red][TRGB][bottom])
	# print(min_t)
	# print(f"minimum exposure time for SNR = 10: {min_t:.2e}s ({min_t/3600:.2e}hrs)")
	
	SNR = {}

	for f in filters:


		SNR[f] = np.ones(len(m_app[red][TRGB]))*10

	"""					

	# this needs some improvement
	pklLoc = f'/Users/harukayoshino/Documents/uvic/masters1/pickles/SNR_{min_SNR}'
	if not os.path.exists(pklLoc): os.makedirs(pklLoc)

	pickleFile = pklLoc + f"/pls_use_this_{bands}.pkl"



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
	"""




	# plot snr vs mag for ALL exposure times for each band (g and i for now)
	# send to pat asap
	# look into roman (how to generate observation) make snr vs mag plot by hand 
	# exposure time needed to reach SNR
	# roman F062 is close to R band? but not exactly matching
	# roman observations in optical or IR?

	##########################
	# APPLY IMF TO ISOCHRONE #
	##########################


	# OBS_SP_RED, OBS_SP_GREEN, SP_FINAL_MASS 
	sp_dictionary = make_sp(red, green, m_app[red][TRGB], 
		m_app[green][TRGB], 1/SNR[red], 1/SNR[green], 
		isochrone_dictionary["final_mass"][TRGB], N_STARS[TRGB])

	if plot:
		plt.scatter(sp_dictionary[green] - sp_dictionary[red], sp_dictionary[red], s=1)
		plt.gca().invert_yaxis()
		plt.show()
		plt.close()


	N_TOT = np.sum(N_STARS[TRGB]) #len(sp_dictionary[red])
	print(N_TOT)
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


	print(tabulate(mag_comp, headers=mag_comp_headers, tablefmt="grid", disable_numparse=True))



	##################
	# SERSIC PROFILE #
	##################
	
	# r_input, I_input_fit, r_obs, I_obs, mu_obs
	luminosity_dictionary = get_luminosity_profile(R_ba, ba, n, Re, b_n, R_MAX, D, sp_dictionary[green])
	

	if plot:
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

		plt.show()
		plt.close()

		# exit()




	outloc = f"/Users/harukayoshino/Documents/uvic/masters1/simulations/{obs_telescope}/maglim_{mag_limit}/{bands}/{green}_{red}/"
	if not os.path.exists(outloc): os.makedirs(outloc)

	with open(outloc+f"{ID}_dictionaries.pkl", "wb") as f:
			pickle.dump([D, AGE, FEH, obs_telescope, red, green, mag_limit, #pickleFile, 
				galaxy_parameters, sp_dictionary, pos_dictionary, luminosity_dictionary], f)

###
# SAVE ALL THE DICTIONARIES BY DUMPING THEM IN A PICKLE FILE OR SOMETHING
# THINK OF A GOOD NAME FOR IT - NEED TO MAKE IT MY UNIQUE CATALOGUE NAME
# VCC929-TRY1, VCC929-TRY2 ...
# MAKE A CSV FILE THAT GIVES THE CATALOGUE NAME AND IDENTIFYING FEATURES (DISTANCE, AGE, FE/H OF ISOCHRONE, EXPOSURE TIME, ETC.)
###


	if save_pdf:
		plotter(ID, galaxy_parameters["RA"], galaxy_parameters["DEC"], 
				galaxy_parameters["Re_arcsec"], n, PA_rad, Re, FEH, AGE, D, SNR, 
				filters, red, green,
				m_app, isochrone_dictionary, sp_dictionary, TRGB, N_STARS, dM, 
				pos_dictionary, luminosity_dictionary, outloc+f"{ID}_sim.pdf")



# for i in range(1,40):
main(plot=False, save_pdf = True, ID_num=1362, lookup = False, bands = "IR", user_min_t = [0,0,0,0,0])

