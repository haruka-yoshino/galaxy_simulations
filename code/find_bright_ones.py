import numpy as np
from astropy.table import Table
from astropy.modeling.models import Sersic1D
import matplotlib.pyplot as plt
import os 
from my_functions import *
from tabulate import tabulate
filename = "/Users/harukayoshino/Documents/uvic/masters1/ngvs catalogue/Galaxy Catalogue/NGVS_CATALOGUE_FINAL_28Mar2026.csv"


tbl = Table.read(filename, format="ascii")


names = np.array(tbl["Nickname"])

col_u = np.array(tbl["u_mag_total"])

col_g = np.array(tbl["g_mag_total"])

col_r = np.array(tbl["r_mag_total"])

col_i = np.array(tbl["i_mag_total"])

col_z = np.array(tbl["z_mag_total"])

Re_g = np.array(tbl["g_re_total"])
mu_e_g = np.array(tbl["g_mue_outer"])
sersic_n = np.array(tbl["g_n_outer"])
PA_g = np.array(tbl["g_av_pa"])
ell = np.array(tbl["g_av_ell"])

GALFIT = np.array(tbl["GALFIT_fit_quality"])
TYPHON = np.array(tbl["TYPHON_nuc_quality"])

detections = (col_g != -100) & (Re_g != -100) & (mu_e_g != -100) & (sersic_n != -100) & (PA_g != -100) & (ell != -100) 

mask_g = (col_g < 15) & detections # & (GALFIT == "good")

print(len(col_g))

print(len(col_g[mask_g]))

# print(names[mask_g])

# exit()

# print(names[mask_g])

obs_telescope = "HST_WFC3"
bands = "ugriz"
D = 16.5*10**6
AGE = 10.3
FEH = -2.0
for ID in names[mask_g]:

	print(ID)

	outloc = f"/Users/harukayoshino/Documents/uvic/masters1/simulations/{obs_telescope}/{bands}/{green}_{red}/{ID}"

	if not os.path.exists(outloc): os.makedirs(outloc)

	galaxy_parameters = get_parameters(filename, ID)

	ba = 1 - galaxy_parameters["ell"]

	Re = D*np.tan(galaxy_parameters["Re_arcsec"]*np.pi/3600/180)

	data = [["ID", ID, "-"], 
	["Distance", f"{D/(10**6):.1f}", "Mpc"], 
	["Age", f"{10**AGE / (10**9):.1f}", "Gyr"], 
	["[Fe/H]", FEH, "-"], 
	["RA", galaxy_parameters["RA"], "deg"], ["DEC", galaxy_parameters["DEC"], "deg"], 
	["GALFIT", galaxy_parameters["GALFIT"], "-"], ["TYPHON", galaxy_parameters["TYPHON"], "-"], 
	["g", galaxy_parameters["mag_total"]["CFHT_g"], "mag"], 
	["R_e", galaxy_parameters["Re_arcsec"], "arcsec"], ["-", Re, "pc"], 
	["mu_{e,g}", galaxy_parameters["mu_eg"], "mag/arcsec^2"], 
	["n", galaxy_parameters["sersic_n"], "-"], 
	["PA", galaxy_parameters["PA"], "deg"], ["-", galaxy_parameters["PA"]*np.pi/180, "rad"], 
	["ellipticity", galaxy_parameters["ell"], "-"], 
	["b/a", ba, "-"]]

	headers = ["Parameter", "Value", "Units"]

	with open(outloc+f"/{ID}_params.txt", "w") as f:
	    f.write(tabulate(data, headers=headers, tablefmt="grid"))

