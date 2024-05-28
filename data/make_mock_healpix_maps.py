import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import healpy as hp
import matplotlib.pyplot as plt
#from healpy.newvisufunc import projview, newprojplot
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

zmin = 0.60
zmax = 0.70


hemi = 'North'
type = 'galaxy'

nside = 128

# Load data and convert into space-separated x, y, z, w in fiducial cosmology
map_dat = np.zeros(12*nside**2)
for i in range(2048):
	if i%32 == rank:
		#data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data
		#data = data[(data['z'] >= zmin) & (data['z'] <= zmax)]
		#data_wt = data['weight_fkp'] * data['weight_systot'] * (data['weight_cp'] + data['weight_noz']- 1)
		data = np.loadtxt('patchy_mocks/Patchy-Mocks-DR12NGC-COMPSAM_V6C_%04d.dat' % (i + 1))
		data = data[(data[:,2] >= zmin) & (data[:,2] <= zmax)]
		data_wt = data[:,6] * data[:,7]/(1 + 10000*data[:,4])

		pixd = hp.ang2pix(nside, data[:,0], data[:,1], lonlat=True)
		map_dat_ind = np.bincount(pixd, minlength=12*nside**2, weights=data_wt)
		hp.write_map('patchy_mocks/healpix_maps/Patchy-Mocks-DR12NGC-COMPSAM_V6C_%04d_%.2f_%.2f_nside128.fits' % (i+1, zmin, zmax), map_dat_ind, overwrite=True)
		map_dat += map_dat_ind
		print(i)




'''#type = 'random0'

#random = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data
#random = random[(random['z'] >= zmin) & (random['z'] <= zmax)]
#random_wt = random['weight_fkp']
random = np.loadtxt('patchy_mocks/Patchy-Mocks-Randoms-DR12NGC-COMPSAM_V6C_x100.dat')
random = random[(random[:,2] >= zmin) & (random[:,2] <= zmax)]
random_wt = random[:,5] * random[:,6]/(1 + 10000*random[:,3])

pixr = hp.ang2pix(nside, random[:,0],random[:,1],lonlat=True)
map_ran = np.bincount(pixr, minlength=12*nside**2,weights=random_wt)


#hp.visufunc.mollview(map,rot=(90,0))
#hp.graticule()
d_minus_r = map_dat-np.sum(map_dat)/np.sum(map_ran) * map_ran
#map_ran_NGC = np.copy(map_ran)


d_minus_r[map_ran == 0] = np.nan
projview(
    d_minus_r, graticule=True, graticule_labels=True, projection_type="mollweide", longitude_grid_spacing = 20, latitude_grid_spacing=10 ,rot=(90,0), xsize=2000
)

plt.savefig('D-R_Patchy-Mocks-DR12NGC-COMPSAM_V6C_zcut_%.2f_%.2f.png' % (zmin, zmax))'''
