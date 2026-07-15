import glob
import numpy as np 

import matplotlib.pyplot as plt 

from scipy.special import gamma, gammainc

import re

from astropy.table import Table
from scipy.interpolate import interp1d



def get_parameters(_file, _ID):
	"""
	get_parameters(_file, _ID)

	This function gets the parameter values from the NGVS catalogue.

	_file (str): path to the catalogue file

	_ID (str): "Nickname" in the NGVS catalogue

	returns RA, DEC, total magnitude in g band, effective radius in g band, 
		surface brightness at effective radius in g band, Sersic n, position angle, 
		and ellipticity

	"""
	_galaxy_parameters = {}

	_tbl = Table.read(_file, format="ascii")



	_row = _tbl[_tbl["Nickname"] == _ID]

	_param_names = ['NGVS_name', 'NGVS_ra(deg)', 'NGVS_dec(deg)', 
	'GALFIT_fit_quality', 'TYPHON_nuc_quality',
	'u_mag_total', 'g_mag_total', 'r_mag_total', 'i_mag_total', 'z_mag_total', 
	'g_re_total', 'g_mue_outer', 'g_n_outer', 'g_av_pa', 'g_av_ell']

	

	NGVS_NAME, RA, DEC, galfit_fit, typhon_fit, u_mag_tot, g_mag_tot, r_mag_tot, i_mag_tot, z_mag_tot, Re_arcsec, mu_eg, n, PA, ell = _row[_param_names][0]

	_galaxy_parameters = {"NGVS_NAME": NGVS_NAME,
	"RA": RA,
	"DEC": DEC,
	"GALFIT": galfit_fit,
	"TYPHON": typhon_fit,
	"mag_total": {"CFHT_u":u_mag_tot, 
				"CFHT_g": g_mag_tot, 
				"CFHT_r": r_mag_tot, 
				"CFHT_i_new": i_mag_tot, 
				"CFHT_z": z_mag_tot},
	"Re_arcsec": Re_arcsec,
	"mu_eg": mu_eg,
	"sersic_n": n,
	"PA": PA,
	"ell": ell
	}
	return _galaxy_parameters



def read_iso_file(_iso_filename):
	"""
	read_iso_file(_iso_filename)

	This function gets the MIST isochrones for a specified [Fe/H].

	Pretty sure this is from the MIST website.

	_iso_filename (str): path to the isochrone file

	returns MIST-MESA version, elemental abundances, stellar rotation, list of ages,
		number of ages, header (with variable names), and data (with variable values)
	"""

	#open file and read it in
	with open(_iso_filename) as _f:
		_content = [_line.split() for _line in _f]
	_version = {'MIST': _content[0][-1], 'MESA': _content[1][-1]}

	_abun = {_content[4][_i]:float(_content[5][_i]) for _i in range(1,5)}

	_rot = float(_content[5][-1])

	_num_ages = int(_content[7][-1])

	#read one block for each isochrone
	_iso_set = []
	_ages = []
	_counter = 10

	for _i_age in range(_num_ages):
		#grab info for each isochrone
		_num_eeps = int(_content[_counter][-2])
		_num_cols = int(_content[_counter][-1])
		_hdr_list = _content[_counter+2][1:]
		_formats = tuple([np.int32]+[np.float64 for _i in range(_num_cols-1)])
		_iso = np.zeros((_num_eeps),{'names':tuple(_hdr_list),'formats':tuple(_formats)})
		#read through EEPs for each isochrone
		for _eep in range(_num_eeps):
			_iso_chunk = _content[3+_counter+_eep]
			_iso[_eep]=tuple(_iso_chunk)
		_iso_set.append(_iso)
		_ages.append(_iso[0][1])
		_counter+= 3+_num_eeps+2

	_ages = np.array([round(_v, 2) for _v in _ages])
	return _version, _abun, _rot, _ages, _num_ages, _hdr_list, _iso_set  


def mag_conv(_telescope, _conversion_dictionary, _f):

	if _telescope == "WFIRST":
		_f = "WFIRST_"+_f

	if _conversion_dictionary[_f]["system"] == "Vega":
		return _conversion_dictionary[_f]["mag(Vega/AB)"]

	elif _conversion_dictionary[_f]["system"] == 'AB':
		return 0


def name_change(_telescope, _spec_band, _v):
	if _telescope == "HST_WFC3":
		if _spec_band == "IR":
			v_i = f"WFC3_IR_{_v}"
		elif _spec_band == "ugriz":
			v_i = f"WFC3_UVIS_{_v}"
	elif _telescope == "HST_ACS":
		v_i = f"ACS_WFC_{_v}"
	else:
		v_i = _v

	return v_i


def get_isochrone(_age=-99.0, _feh=-99.0,  #_isochrone_dictionary = None, 
	_telescope=None, _spec_band = None, _filters=None, plot=False, _xlabel=None, _ylabel=None):
	"""
	get_isochrone(_age=-99.0, _feh=-99.0, _telescope=None, _filters=None)

	This function finds an isochrone that matches the user-defined age and metallicity.

	If the parameters (_age and _feh) are not defined before the function is called (left 
	as -99.0), it will ask for user input.

	If the parameters are defined and called into the function, it will skip the step 
	requiring user input.
	"""
	mag_conversion_file = "/Users/harukayoshino/Documents/uvic/masters1/isochrones/mag_conversion.txt"
	
	tbl = Table.read(mag_conversion_file, format="ascii")
	
	# First column becomes top-level keys
	first_col = tbl.colnames[0]
	other_cols = tbl.colnames[1:]

	_conversion_dictionary = {
		row[first_col]: {col: row[col] for col in other_cols}
		for row in tbl
	}

	_path = f"/Users/harukayoshino/Documents/uvic/masters1/isochrones/*_{_telescope}/*"
	_temp = glob.glob(_path)

	_feh_list=[]

	for _t in _temp:
		_match = re.search(r'feh_([mp]?\d+\.\d+)', _t)
		_val = _match.group(1)
		_val = float(_val.replace('m', '-').replace('p', ''))
		_feh_list.append(_val)

	_feh_list = np.array(_feh_list)

	while _feh == -99.0:
		print("[Fe/H] options:", np.array(sorted(_feh_list)))

		_feh_try = float(input("[Fe/H]: "))

		if np.isin(_feh_try, _feh_list):
			_feh = _feh_try

		else:
			input(f"[Fe/H] = {_feh_try} not available")

	_indx = np.where(_feh_list == _feh)[0][0]

	_iso_temp = _temp[_indx]

	_version, _abun, _rot, _age_list, _num_ages, _hdr_list, _iso_set  = read_iso_file(_iso_temp)


	while _age == -99.0:

		print("log10(Ages) options:", _age_list)

		_age_try = float(input("log10(Age): "))

		if np.isin(_age_try, _age_list):
			_age = _age_try

		else:
			input(f"log10(Age) = {_age_try} not available")

	_indx = np.where(_age_list == _age)[0][0]

	one_iso = _iso_set[_indx]


	mapping = {"init_mass": "initial_mass",
	"final_mass": "star_mass",
	"M_abs": _filters,
	"T_eff": "log_Teff",
	"star_type": "phase"
	}

	# init_mass_indx = _hdr_list.index('initial_mass')
	# final_mass_indx = _hdr_list.index('star_mass')

	_isochrone_dictionary = {} #"AGE": _age, "FEH": _feh
	for label, var in mapping.items():
		if label == "M_abs":
			mags_dict = {}

			for v in mapping[label]:
				v_i = name_change(_telescope, _spec_band, v)

				m_Vega_AB = mag_conv(_telescope, _conversion_dictionary, v_i)

				_var_indx = _hdr_list.index(v_i)
				mags_dict[v] = np.array([t[_var_indx] + m_Vega_AB for t in one_iso]) 


			_isochrone_dictionary[label] = mags_dict# = np.array([t[_var_indx] for t in one_iso]) 
	
		else:
			_var_indx = _hdr_list.index(var)
			_isochrone_dictionary[label] = np.array([t[_var_indx] for t in one_iso]) 
	
	_isochrone_dictionary["T_eff"] = 10**_isochrone_dictionary["T_eff"]
	_isochrone_dictionary["RGB"] = (_isochrone_dictionary["star_type"] == 2.0)

	# if plot:
	# 	print("insert updated plots here")


	return _isochrone_dictionary 


def get_CFHT_mags(_age=-99.0, _feh=-99.0,  #_isochrone_dictionary = None, 
	_telescope=None, _spec_band = None, _filters=None, plot=False, _xlabel=None, _ylabel=None):
	"""
	get_isochrone(_age=-99.0, _feh=-99.0, _telescope=None, _filters=None)

	This function finds an isochrone that matches the user-defined age and metallicity.

	If the parameters (_age and _feh) are not defined before the function is called (left 
	as -99.0), it will ask for user input.

	If the parameters are defined and called into the function, it will skip the step 
	requiring user input.
	"""
	mag_conversion_file = "/Users/harukayoshino/Documents/uvic/masters1/isochrones/mag_conversion.txt"
	
	tbl = Table.read(mag_conversion_file, format="ascii")
	
	# First column becomes top-level keys
	first_col = tbl.colnames[0]
	other_cols = tbl.colnames[1:]

	_conversion_dictionary = {
		row[first_col]: {col: row[col] for col in other_cols}
		for row in tbl
	}

	_path = f"/Users/harukayoshino/Documents/uvic/masters1/isochrones/*_{_telescope}/*"
	_temp = glob.glob(_path)

	_feh_list=[]

	for _t in _temp:
		_match = re.search(r'feh_([mp]?\d+\.\d+)', _t)
		_val = _match.group(1)
		_val = float(_val.replace('m', '-').replace('p', ''))
		_feh_list.append(_val)

	_feh_list = np.array(_feh_list)

	while _feh == -99.0:
		print("[Fe/H] options:", np.array(sorted(_feh_list)))

		_feh_try = float(input("[Fe/H]: "))

		if np.isin(_feh_try, _feh_list):
			_feh = _feh_try

		else:
			input(f"[Fe/H] = {_feh_try} not available")

	_indx = np.where(_feh_list == _feh)[0][0]

	_iso_temp = _temp[_indx]

	_version, _abun, _rot, _age_list, _num_ages, _hdr_list, _iso_set  = read_iso_file(_iso_temp)


	while _age == -99.0:

		print("log10(Ages) options:", _age_list)

		_age_try = float(input("log10(Age): "))

		if np.isin(_age_try, _age_list):
			_age = _age_try

		else:
			input(f"log10(Age) = {_age_try} not available")

	_indx = np.where(_age_list == _age)[0][0]

	one_iso = _iso_set[_indx]

	mags_dict = {}

	for v in _filters:
		v_i = name_change(_telescope, _spec_band, v)

		m_Vega_AB = mag_conv(_telescope, _conversion_dictionary, v_i)

		_var_indx = _hdr_list.index(v_i)
		mags_dict[v] = np.array([t[_var_indx] + m_Vega_AB for t in one_iso]) 


	return mags_dict 


def extract_TRGB(_m_green, _m_red, mag_limit = 28.0, color_range=[0.0,3.0]):
    
    # observational magnitude limit 
    # we can't observe anything below this magnitude with the telescope used
    _mag_lower = np.where(_m_red <= mag_limit)[0]#np.min(_m_red)+2)[0]

    # color index range to exclude unwanted stars
    _col_upper = np.where((_m_green - _m_red) <= color_range[1])[0]
    _col_lower = np.where((_m_green - _m_red) >= color_range[0])[0]

	# put it all together
    _mag_indx = list(set(_mag_lower) & set(_col_upper) & set(_col_lower)) 
    
    return _mag_indx 
    

def salpeter_imf(_init_mass):
	"""
	Salpeter IMF
	"""
	return _init_mass**(-2.35)


def flux(_mag):
	return 10**(-0.4*(_mag))

def luminosity(_mag, _d):
    return 4 * np.pi * _d**2 * flux(_mag)


def magnitude(_flux):
	return -2.5*np.log10(_flux)


def norm_imf(_g, _init_mass, _m_green, _D, _dM, _imf):
    Lg_tot = luminosity(_g, _D) 
    Lg_M = luminosity(_m_green, _D) 
    
    _norm = Lg_tot / np.sum(Lg_M*_imf(_init_mass)*_dM)

    return _norm


def make_sp(_red, _green, _m_red, _m_green, _m_red_err, _m_green_err, _final_mass, _n_stars):

	_sp_red = np.repeat(_m_red, _n_stars)
	_sp_green = np.repeat(_m_green, _n_stars)

	_sp_red_err = np.repeat(_m_red_err, _n_stars)
	_sp_green_err = np.repeat(_m_green_err, _n_stars)

	_sp_final_mass = np.repeat(_final_mass, _n_stars)


	_obs_sp_red = np.random.normal(loc=_sp_red, scale=_sp_red_err) 
	_obs_sp_green = np.random.normal(loc=_sp_green, scale=_sp_green_err) 


	_sp_dictionary = {_red: _obs_sp_red, _green: _obs_sp_green, 
		"final_mass": _sp_final_mass
		}

	return _sp_dictionary #_sp_red, _sp_green, _sp_red_err, _sp_green_err, _sp_final_mass



def sample_sersic(_n, _Re, _N, _bn, r_max):

	_Rlist = np.linspace(0, np.log10(r_max), 10000)
	_I = sersic_profile(_n, 10**_Rlist, _Re, _bn)
	_P = _I * (10**_Rlist)**2
	_P /= np.sum(_P)
	sampled_R = np.random.choice(10**_Rlist, size=_N, p=_P)

	return sampled_R


def get_cartesian(_ba, _R, _theta):

    _x = _ba*_R * np.sin(_theta)
    _y = _R * np.cos(_theta)
    
    return _x, _y


def rotate_clockwise(_x, _y, _PA):
	_xrot = _x*np.cos(_PA) + _y*np.sin(_PA)
	_yrot = -_x*np.sin(_PA) + _y*np.cos(_PA)

	return _xrot, _yrot

def get_positions(_n_tot, _ba, _R, _PA_rad):

	_theta = np.random.uniform(0, 2*np.pi, _n_tot)

	_xba, _yba = get_cartesian(_ba, _R, _theta)

	_xrot, _yrot = rotate_clockwise(_xba, _yba, _PA_rad)

	_R_ba = np.sqrt((_xba)**2 + (_yba)**2)


	return _xrot, _yrot, _R_ba



def get_ellipses(_R_circ, _ba, _PA_rad):

	_x_circ = []
	_y_circ = []

	_theta_circ = np.linspace(0, 2*np.pi, 100) 

	for r in _R_circ:
		_xba, _yba = get_cartesian(_ba, r, _theta_circ)

		_xrot, _yrot = rotate_clockwise(_xba, _yba, _PA_rad)


		_x_circ.append(_xrot)
		_y_circ.append(_yrot)

	return _x_circ, _y_circ


def sersic_profile(_n, _R, _Re, _bn):
    # _bn = 2*_n - 1/3 + 4/(405*_n) + 46/(25515*_n**2) + 131/(1148175*_n**3) - 2194697/(30690717750*_n**4)

    _I = np.exp(-_bn*((_R / _Re)**(1/_n) - 1))

    return _I


def get_annuli(_edges, _ba):

	# get area of ellipse at each edge
	_ellipses = np.pi*_ba*_edges**2

	# get area of each annulus
	_annuli = _ellipses[1:] - _ellipses[:-1]

	return _annuli



def get_luminosity_profile(_R_ba, _ba, _sersic_n, _Re, _b_n, r_max, _D, sp_mag_g):
	"""
	this function needs more construction work pls
	"""


	# bin radii
	counts, edges_pc = np.histogram(_R_ba, bins=10**np.linspace(0, np.log10(r_max), 20)) # parsecs

	# get area of annuli in pc^2
	annuli_pc = get_annuli(edges_pc, _ba)


	# find the bin in which each star belongs to
	bin_indices = np.digitize(_R_ba, edges_pc) - 1
	bin_indices = np.clip(bin_indices, 0, len(counts)-1)



	# find the average radius of the stars in each bin
	centers_avg = np.array([np.mean(_R_ba[bin_indices == i]) for i in range(len(counts))])


	# find the total g band magnitude per bin
	mag_g_total = np.array([magnitude(np.sum(flux(sp_mag_g[bin_indices == i]))) for i in range(len(counts))])



	# get area of annuli in arcsec^2
	annuli_as = get_annuli(pc_to_arcsec(edges_pc, _D), _ba)


	# get magnitude / arcsec^2
	mu_g_binned = mag_g_total + 2.5 * np.log10(annuli_as)



	# get binned radius and surface density in each bin
	r_binned = np.log10(centers_avg)
	I_binned = np.log10(counts / annuli_pc)


	# mask out infinite values
	mask_binned = np.isfinite(r_binned) & np.isfinite(I_binned) & np.isfinite(mu_g_binned)


	# line of best fit part

	r_bf = np.linspace(0, np.log10(r_max), 1000) #np.log10(np.sort(_R_ba))

	I_bf = np.log10(sersic_profile(_sersic_n, 10**r_bf, _Re, _b_n))


	# mask out infinite values
	mask_bf = np.isfinite(r_bf) & np.isfinite(I_bf)




	# interpolate between all the sampled radii
	interp_model = interp1d(r_bf[mask_bf], I_bf[mask_bf], kind='linear', fill_value="extrapolate")

	# get where the line of best fit should be
	I_input_binned = interp_model(r_binned[mask_binned]) # what it should actually be



	shift_bf = np.mean(I_binned[mask_binned] - I_input_binned) # get the mean vertical shift from what it should actually be

	I_bf_fit = I_bf[mask_bf] + shift_bf



	_luminosity_dictionary = { "r_bf": r_bf[mask_bf], "I_bf": I_bf_fit[mask_bf], 
				"r_binned": r_binned[mask_binned], "I_binned": I_binned[mask_binned], 
				"mu_e": mu_g_binned[mask_binned] }

	return _luminosity_dictionary #r_bf[mask_bf], I_bf_fit[mask_bf], r_binned[mask_binned], I_binned[mask_binned], mu_g_binned[mask_binned]




# conversion functions
def pc_to_arcsec(_r, _D):
    return np.arctan(_r / _D) * 3600 * 180 / np.pi

def arcsec_to_pc(_theta, _D):
    return _D * np.tan(_theta * np.pi / 3600 / 180)



def get_detections(_filename, magrange):# max_mag = False, min_mag = False):
	min_mag, max_mag = magrange
	print(min_mag, max_mag)
	_tbl = Table.read(_filename, format="ascii")


	_names = np.array(_tbl["Nickname"])

	_col_g = np.array(_tbl["g_mag_total"])

	_col_Re_g = np.array(_tbl["g_re_total"])
	_col_sersic_n = np.array(_tbl["g_n_outer"])
	_col_PA_g = np.array(_tbl["g_av_pa"])
	_col_ell = np.array(_tbl["g_av_ell"])


	_detections = (_col_g != -100) & (_col_Re_g != -100) & (_col_sersic_n != -100) & (_col_PA_g != -100) & (_col_ell != -100) 

	_mask_g = _detections 

	if max_mag:
		_mask_g = _mask_g & (_col_g < max_mag) 
	

	if min_mag:
		_mask_g = _mask_g & (_col_g >= min_mag) 


	return _names[_mask_g]#, _col_Re_g[_mask_g], _col_PA_g[_mask_g], 1-_col_ell[_mask_g]

