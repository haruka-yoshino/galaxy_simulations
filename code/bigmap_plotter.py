import matplotlib.pyplot as plt
import mpl_scatter_density
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from my_functions import *
from astropy.io import ascii

import numpy as np 

import os 
import pickle
from astropy.table import Table
from astropy.io import fits


plt.style.use("/Users/harukayoshino/Documents/uvic/masters1/code/paper_plots.mplstyle")


# obs_telescope = "HST_WFC3"
obs_telescope = "WFIRST"
bands = "IR" #"ugriz"
mag_limit =  27.52 #27.9 #

red = "H158" #"F140W" ##"CFHT_i_new" "F160W" #
green = "J129" #"F105W" #"CFHT_g" "F125W" #

maps_loc = f"/Users/harukayoshino/Documents/uvic/masters1/simulations/{obs_telescope}/maglim_{mag_limit}/{bands}/{green}_{red}/exp_time_10hr/NGVS_full_maps"
if not os.path.exists(maps_loc): os.makedirs(maps_loc)


filename = "/Users/harukayoshino/Documents/uvic/masters1/ngvs catalogue/Galaxy Catalogue/NGVS_CATALOGUE_FINAL_28Mar2026.csv"

footprint = "/Users/harukayoshino/Documents/uvic/masters1/ngvs catalogue/footprint_dictionary.pkl"


with open(footprint, 'rb') as f:  # Python 3
	footprint_dictionary = pickle.load(f)

names = get_detections(filename, [8,22])

# # = ["NGVS1926", "NGVS2115"]
# m49 = np.where(names == "NGVS1926")[0][0]
# m87 = np.where(names == "NGVS2115")[0][0]


xmin = 181 #180.93775 #- 2
xmax = 194 #193.3790833 #+ 2
ymin = 4 #4.8495833 #- 2
ymax = 19 #27.6172222 #+ 2

x_range = xmax - xmin
y_range = ymax - ymin




# print(col_RA[m49], col_DEC[m49])
# exit()

# exit()



all_RA = []
all_DEC = []

all_red = []
all_green = []

center_RA = []
center_DEC = []

# one = 1
# print(names[mask_g][one])
# NSTARS = 0


_theta_circ = np.linspace(0, 2*np.pi, 100) 


XROT_circ = []
YROT_circ = []


max_R = []
# all_Re = []

# all_Re_arcsec = []
# all_Re_pc = []
# all_ba = []
# all_PA = []

for ID in names:#[[m49, m87]]: #[names[mask_g][one]]: #i in range(0, 201): #

	

	outloc = f"/Users/harukayoshino/Documents/uvic/masters1/simulations/{obs_telescope}/maglim_{mag_limit}/{bands}/{green}_{red}/exp_time_10hr/{ID}"
	if os.path.exists(outloc+f"/{ID}_dictionaries.pkl"):
		# if it pre-exists, open it
		with open(outloc+f"/{ID}_dictionaries.pkl", 'rb') as f:  # Python 3	
			D, AGE, FEH, obs_telescope, red, green, mag_limit, pickleFile, galaxy_parameters, sp_dictionary, pos_dictionary, luminosity_dictionary = pickle.load(f)

		mags_red = sp_dictionary[red]
		mags_green = sp_dictionary[green]

	else: 
		continue


	if os.path.exists(outloc+f"/{ID}_locs.pkl"):
		# if it pre-exists, open it
		with open(outloc+f"/{ID}_locs.pkl", 'rb') as f:  # Python 3
			ID, _RA0, _DEC0, _Re_arcsec, _xrot, _yrot, _Re_pc, _ba, _PA = pickle.load(f)


		
		RA_pos = np.arctan(_xrot/(16.5*10**6)) * 180 / np.pi +_RA0
		DEC_pos = np.arctan(_yrot/(16.5*10**6)) * 180 / np.pi +_DEC0

		max_R.append(np.max(np.sqrt((_xrot)**2 + (_yrot)**2)))


		mask = (RA_pos >= xmin) & (RA_pos <= xmax) & (DEC_pos >= ymin) & (DEC_pos <= ymax)



		all_RA.append(RA_pos[mask])
		all_DEC.append(DEC_pos[mask])
		all_red.append(mags_red[mask])
		all_green.append(mags_green[mask])


		# all_Re_arcsec.append(_Re_arcsec)
		# all_Re_pc.append(_Re_pc)
		# all_ba.append(_ba)
		# all_PA.append(_PA)

		center_RA.append(_RA0)
		center_DEC.append(_DEC0)




		# _Re = all_Re_arcsec / 3600 # all_Re[i]
		# _ba = all_ba
		# _PA = PAs[j] * np.pi / 180




		# _RA0 = center_RA[i]
		# _DEC0 = center_DEC[i]



		# _xell, _yell = get_ellipses([3*_Re_arcsec], _ba, _PA)

		# XROT_circ.append(_xell+_RA0)
		# YROT_circ.append(_yell+_DEC0)

	else:
		continue




# for i, n in enumerate(names):#[[m49, m87]]):

# 	_Re = all_Re_arcsec / 3600 # all_Re[i]
# 	_ba = all_ba
# 	_PA = PAs[j] * np.pi / 180

# 	_RA0 = center_RA[i]
# 	_DEC0 = center_DEC[i]

# 	_xell, _yell = get_ellipses([3*_Re], _ba, _PA)

# 	XROT_circ.append(_xell+_RA0)
# 	YROT_circ.append(_yell+_DEC0)





all_RA = np.concatenate(all_RA)
all_DEC = np.concatenate(all_DEC)
all_red = np.concatenate(all_red)
all_green = np.concatenate(all_green)

print(f"{len(all_RA):.2e}")




# xmin = np.min(all_RA) # 180.93775 #- 2
# xmax = np.max(all_RA) # 193.3790833 #+ 2
# ymin = np.min(all_DEC) # 4.8495833 #- 2
# ymax = np.max(all_DEC) # 18.5 # 27.6172222 #+ 2


# _virgo_center_RA = np.median(all_RA)
# _virgo_center_DEC = np.median(all_DEC)

# _Redge = np.arctan((1.5*10**6) / (16.5*10**6)) * 180 / np.pi

# print(_Redge)

# _xedge = _Redge*np.sin(_theta_circ) + _virgo_center_RA
# _yedge = _Redge*np.cos(_theta_circ) + _virgo_center_DEC


bin_width = 0.001 #0005 #10**2 # pc resolution 
# _square = np.max([-np.min(_xrot), np.max(_xrot), -np.min(_yrot), np.max(_yrot)]) / (10**3)

bins_x = np.arange(xmin, xmax + bin_width, bin_width)
bins_y = np.arange(ymin, ymax + bin_width, bin_width)

print(len(bins_x), len(bins_y))

H, xedges, yedges = np.histogram2d(all_RA, all_DEC, bins = [bins_x, bins_y], weights=flux(all_red))


# dx = xedges[1] - xedges[0]
# dy = yedges[1] - yedges[0]
# density = H.T / (dx * dy)

# xcenters = (xedges[:-1] + xedges[1:]) / 2
# ycenters = (yedges[:-1] + yedges[1:]) / 2

# fits.writeto(f"/Users/harukayoshino/Documents/uvic/masters1/simulations/{obs_telescope}/{bands}/{green}_{red}/NGVS_full_maps/fits_image_try_zoomed_highres.fits", density, overwrite=True)
fig = plt.figure(figsize=(12+2, 12 * y_range / x_range))#figsize = (ymax - ymin, xmax - xmin))

ax = fig.add_subplot(projection='scatter_density')


# cmap = plt.get_cmap('plasma')
# darkest_color = cmap(0.0) 

 
# ax.set_facecolor(darkest_color)

# plt.close()
im = ax.imshow(magnitude(H.T), origin='lower', 
	cmap = "plasma_r",  aspect="equal", extent = [xmin, xmax, ymin, ymax],
	vmin = 14.5, vmax = 28.5)




# density22 = ax.scatter_density(all_RA, all_DEC, cmap="plasma")

# fig.canvas.draw()
# vals22 = density22.get_array()

# print(np.shape(vals22))
# vals22 = vals22[np.isfinite(vals22) & (vals22 > 0)]

# vmin22 = 1.0 #np.percentile(vals22, 0.01)
# vmax22 = 700728.40 #np.percentile(vals22, 99.99)
# print(vmin22, vmax22)
# density22.set_norm(LogNorm(vmin=vmin22, vmax=vmax22))
cbar = fig.colorbar(im, label="total magnitude in bin")
cbar.ax.invert_yaxis()
# axc = inset_axes(ax, width='5%', height='100%', loc='upper center', bbox_to_anchor=(1.4, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=0)
# axcb = fig.colorbar(density22, cax=axc, label="number / square degree")



# for i in range(len(names)):
# 	if len(XROT_circ[i]) != 0:
# 		ax.scatter(XROT_circ[i], YROT_circ[i], color="green", s=0.05)#, label = ["M49", "M87"][i])
# 	else: 
# 		pass
# # for i in range(0, len(center_RA)):
# # 	ax.scatter(center_RA[i], center_DEC[i], color="red", s=col_g[mask_g][i])
# # 	ax.text(center_RA[i], center_DEC[i], names[mask_g][i], fontsize=3)
# # ax.scatter(center_RA[m49], center_DEC[m49], color="green", s=5)#marker="*", s=20)
# # ax.scatter(center_RA[m87], center_DEC[m87], color="black", s=5)#marker="*", s=20)

# # ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, fancybox=True, shadow=True)
ax.plot(footprint_dictionary["ra"], footprint_dictionary["dec"], color="green")
ax.set_xlabel("RA (deg)")
ax.set_ylabel("DEC (deg)")
ax.set_title(f"Virgo Cluster ({obs_telescope})\n({bin_width*3600} arcsecond resolution)\n")
# ax.axis("equal")

# ax.set_xlim(xmax, xmin)
# ax.set_ylim(ymin, ymax)

ax.invert_xaxis()
plt.savefig(f"{maps_loc}/virgo_resolved_map_{bands}_8-22.pdf",dpi=1200)#_centers_labeled.pdf")
# plt.show()
plt.close()

