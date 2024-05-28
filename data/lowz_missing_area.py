# coding: utf-8
import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import healpy as hp
import matplotlib.pyplot as plt
from healpy.newvisufunc import projview, newprojplot


hemi = 'North'
type = 'galaxy'

# Load data and convert into space-separated x, y, z, w in fiducial cosmology

data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data

nside = 32

pix = hp.ang2pix(nside, data['ra'],data['dec'],lonlat=True)
map = np.bincount(pix, minlength=12*nside**2)
#hp.visufunc.mollview(map,rot=(90,0))
#hp.graticule()
# classic healpy mollweide projections plot with graticule
projview(
    map, min=0, graticule=True, graticule_labels=True, projection_type="mollweide", longitude_grid_spacing = 20, latitude_grid_spacing=10 ,rot=(90,0), xsize=2000
)
plt.savefig('galaxy_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.png' % (hemi))

type = 'random0'

data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data
pix = hp.ang2pix(nside, data['ra'],data['dec'],lonlat=True)
map_ran = np.bincount(pix, minlength=12*nside**2)
#hp.visufunc.mollview(map,rot=(90,0))
#hp.graticule()
projview(
    map_ran, graticule=True, graticule_labels=True, projection_type="mollweide", longitude_grid_spacing = 20, latitude_grid_spacing=10 ,rot=(90,0), xsize=2000
)

plt.savefig('random0_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.png' % (hemi))


projview(
    map/map_ran, min=0.01, max=0.03, graticule=True, graticule_labels=True, projection_type="mollweide", longitude_grid_spacing = 20, latitude_grid_spacing=10 ,rot=(90,0), xsize=2000
)

plt.savefig('galaxy_ran_ratio_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.png' % (hemi))

type = 'galaxy'
data = fits.open('%s_DR12v5_LOWZ_%s.fits.gz' % (type, hemi))[1].data


pix = hp.ang2pix(nside, data['ra'],data['dec'],lonlat=True)
map_lowz = np.bincount(pix, minlength=12*nside**2)
#hp.visufunc.mollview(map,rot=(90,0))
#hp.graticule()
# classic healpy mollweide projections plot with graticule
projview(
    map_lowz, graticule=True, graticule_labels=True, projection_type="mollweide", longitude_grid_spacing = 20, latitude_grid_spacing=10 ,rot=(90,0), xsize=2000
)
plt.savefig('galaxy_DR12v5_LOWZ.png')

print(np.median(map[map_lowz == 0][map[map_lowz==0] > 0]))
print(np.std(map[map_lowz == 0][map[map_lowz==0] > 0]))
print(len(map[map_lowz == 0][map[map_lowz==0] > 0]))


print(np.median(map[map_lowz != 0][map[map_lowz!=0] > 0]))
print(np.std(map[map_lowz != 0][map[map_lowz!=0] > 0]))
print(len(map[map_lowz != 0][map[map_lowz!=0] > 0]))


print(np.median(map[map_lowz == 0][map[map_lowz==0] > 0]/map_ran[map_lowz==0][map[map_lowz==0] > 0]))

print(np.median(map[map_lowz != 0][map[map_lowz!=0] > 0]/map_ran[map_lowz!=0][map[map_lowz!=0] > 0]))



data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data

nside = 16

pix = hp.ang2pix(nside, data['ra'],data['dec'],lonlat=True)
map = np.bincount(pix, minlength=12*nside**2)
#hp.visufunc.mollview(map,rot=(90,0))
#hp.graticule()
# classic healpy mollweide projections plot with graticule
projview(
    map, min=500, graticule=True, graticule_labels=True, projection_type="mollweide", longitude_grid_spacing = 20, latitude_grid_spacing=10 ,rot=(90,0), xsize=2000
)
plt.savefig('galaxy_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70_nside16.png' % (hemi))


data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data

nside = 16

pix = hp.ang2pix(nside, data['ra'][(data['Z'] >= 0.43) & (data['Z'] <= 0.50)],data['dec'][(data['Z'] >= 0.43) & (data['Z'] <= 0.50)],lonlat=True)
map = np.bincount(pix, minlength=12*nside**2)
#hp.visufunc.mollview(map,rot=(90,0))
#hp.graticule()
# classic healpy mollweide projections plot with graticule
projview(
    map, min=500, graticule=True, graticule_labels=True, projection_type="mollweide", longitude_grid_spacing = 20, latitude_grid_spacing=10 ,rot=(90,0), xsize=2000
)
plt.savefig('galaxy_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.50_nside16.png' % (hemi))
data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data

nside = 16

pix = hp.ang2pix(nside, data['ra'][(data['Z'] >= 0.43) & (data['Z'] <= 0.50)],data['dec'][(data['Z'] >= 0.43) & (data['Z'] <= 0.50)],lonlat=True)
map = np.bincount(pix, minlength=12*nside**2)
#hp.visufunc.mollview(map,rot=(90,0))
#hp.graticule()
# classic healpy mollweide projections plot with graticule
projview(
    map, graticule=True, graticule_labels=True, projection_type="mollweide", longitude_grid_spacing = 20, latitude_grid_spacing=10 ,rot=(90,0), xsize=2000
)
plt.savefig('galaxy_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.50_nside16.png' % (hemi))
import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import healpy as hp
import matplotlib.pyplot as plt
from healpy.newvisufunc import projview, newprojplot


hemi = 'North'
type = 'galaxy'

# Load data and convert into space-separated x, y, z, w in fiducial cosmology

data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data

nside = 32

pix = hp.ang2pix(nside, data['ra'],data['dec'],lonlat=True)
map = np.bincount(pix, minlength=12*nside**2)
#hp.visufunc.mollview(map,rot=(90,0))
#hp.graticule()
# classic healpy mollweide projections plot with graticule
projview(
    map, min=0, graticule=True, graticule_labels=True, projection_type="mollweide", longitude_grid_spacing = 20, latitude_grid_spacing=10 ,rot=(90,0), xsize=2000
)
plt.savefig('galaxy_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.png' % (hemi))

type = 'random0'

data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data
pix = hp.ang2pix(nside, data['ra'],data['dec'],lonlat=True)
map_ran = np.bincount(pix, minlength=12*nside**2)
#hp.visufunc.mollview(map,rot=(90,0))
#hp.graticule()
projview(
    map_ran, graticule=True, graticule_labels=True, projection_type="mollweide", longitude_grid_spacing = 20, latitude_grid_spacing=10 ,rot=(90,0), xsize=2000
)

plt.savefig('random0_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.png' % (hemi))


projview(
    map/map_ran, min=0.01, max=0.03, graticule=True, graticule_labels=True, projection_type="mollweide", longitude_grid_spacing = 20, latitude_grid_spacing=10 ,rot=(90,0), xsize=2000
)

plt.savefig('galaxy_ran_ratio_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.png' % (hemi))

type = 'galaxy'
data = fits.open('%s_DR12v5_LOWZ_%s.fits.gz' % (type, hemi))[1].data


pix = hp.ang2pix(nside, data['ra'],data['dec'],lonlat=True)
map_lowz = np.bincount(pix, minlength=12*nside**2)
#hp.visufunc.mollview(map,rot=(90,0))
#hp.graticule()
# classic healpy mollweide projections plot with graticule
projview(
    map_lowz, graticule=True, graticule_labels=True, projection_type="mollweide", longitude_grid_spacing = 20, latitude_grid_spacing=10 ,rot=(90,0), xsize=2000
)
plt.savefig('galaxy_DR12v5_LOWZ.png')

print(np.median(map[map_lowz == 0][map[map_lowz==0] > 0]))
print(np.std(map[map_lowz == 0][map[map_lowz==0] > 0]))
print(len(map[map_lowz == 0][map[map_lowz==0] > 0]))


print(np.median(map[map_lowz != 0][map[map_lowz!=0] > 0]))
print(np.std(map[map_lowz != 0][map[map_lowz!=0] > 0]))
print(len(map[map_lowz != 0][map[map_lowz!=0] > 0]))


print(np.median(map[map_lowz == 0][map[map_lowz==0] > 0]/map_ran[map_lowz==0][map[map_lowz==0] > 0]))

print(np.median(map[map_lowz != 0][map[map_lowz!=0] > 0]/map_ran[map_lowz!=0][map[map_lowz!=0] > 0]))


pixs = np.where((map_lowz ==0) & (map != 0))[0]
pixs
data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data

pix = hp.ang2pix(nside, data['ra'],data['dec'],lonlat=True)
pixs in pix
np.where(pixs in pix)
pixs = np.where((map_lowz ==0) & (map != 0))[0]

data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data

pix = hp.ang2pix(nside, data['ra'],data['dec'],lonlat=True)
zz = []
for i, pp in enumerate(pix):
	if pp in pixs:
		zz.append(data[i]['Z'])
zz = np.array(zz)
np.mean(zz)
np.mean(data['Z'])
np.shape(zz)
np.shape(data['Z'])
mask = np.zeros(12*nside**2)
mask[pixs = 1]
mask[pixs] = 1
hp.write_map('lowz_missing_mask.png',mask)
hp.visufunc.mollview(mask)
plt.savefig('lowz_missing_mask.png')
hp.visufunc.mollview(mask,rot=(90,0))
plt.savefig('lowz_missing_mask.png')
np.shape(np.where(zz < 0.45))
np.shape(np.where(data['Z'] < 0.45))
np.shape(data)
np.shape(zz)
33751./587071.
2774./56078.
plt.ion()
plt.show()
plt.close9"all')
plt.close('all')
plt.figure()
plt.hist(data['Z'],min=0.43,max=0.7,bins=100,label='All data')
plt.hist(data['Z'],range=(0.43,0.7),bins=100,label='All data')
plt.figure()
plt.hist(data['Z'],range=(0.43,0.7),bins=100,label='All data',histtype='step')
plt.hist(zz,range=(0.43,0.7),bins=100,label='LOWZ missing area',histtype='step')
plt.figure()
plt.hist(data['Z'],range=(0.43,0.7),bins=100,label='All data',histtype='step',normed=True)
plt.hist(data['Z'],range=(0.43,0.7),bins=100,label='All data',histtype='step',density=True)
plt.figure()
plt.hist(data['Z'],range=(0.43,0.7),bins=100,label='All data',histtype='step',density=True)
plt.hist(zz,range=(0.43,0.7),bins=100,label='LOWZ missing area',histtype='step',density=True)
plt.legend()
plt.savefig('zdist_lowz_missing_area.png')
