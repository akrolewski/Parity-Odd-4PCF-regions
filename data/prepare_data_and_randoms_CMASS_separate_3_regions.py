import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import time
import sys
import os

#rank = int(sys.argv[1])
#rank = 0

t0 = time.time()

hemi = 'S'

region = 3


# Load data and convert into space-separated x, y, z, w in fiducial cosmology

#data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data
zmin = 0.43
zmax = 0.70

data = fits.open('/gpfs/akrolewski/parity_odd_4pcf/data/galaxy_DR12v5_CMASS_North.fits.gz')[1].data
data = data[(data['Z'] >= zmin) & (data['Z'] <= zmax) ]
if region == 1:
	data = data[(data['Z'] >= zmin) & (data['Z'] <= zmax) 
		& (data['RA']< (155 - 1./12. * data['DEC']))]
elif region == 2:
	data = data[(data['Z'] >= zmin) & (data['Z'] <= zmax) 
		& (data['RA'] > (165 + 1./12. * data['DEC']))
		& (data['RA'] < (205. - 1./12. * data['DEC'])
		)]
elif region == 3:
	data = data[(data['Z'] >= zmin) & (data['Z'] <= zmax) 
		& (data['RA']> (215. + 1./12. * data['DEC'])
		)]

weight = data['WEIGHT_SYSTOT'] * (data['WEIGHT_NOZ'] + data['WEIGHT_CP'] - 1.0) * data['WEIGHT_FKP']
cosmo = FlatLambdaCDM(H0 = 67.6, Om0 = 0.31, m_nu = [0, 0, 0.06] *u.eV, Ob0 = 0.022/(0.676**2.))

rr = cosmo.comoving_distance(data['Z']).value * 0.676

x = rr * np.sin((90. - data['DEC']) * np.pi/180.) * np.cos(data['RA']*np.pi/180.)
y = rr * np.sin((90. - data['DEC']) * np.pi/180.) * np.sin(data['RA']*np.pi/180.)
z = rr * np.cos((90. - data['DEC']) * np.pi/180.)

#if type == 'galaxy':
#weight = data['WEIGHT_SYSTOT'] * (data['WEIGHT_NOZ'] + data['WEIGHT_CP'] - 1.0) * data['WEIGHT_FKP']
sum_wts = np.sum(weight)

out = np.array([x, y, z, weight]).T
np.savetxt('/gpfs/akrolewski/parity_odd_4pcf/data/galaxy_DR12v5_CMASS_North_reg%i_cart.data' % (region),out)
#print(j)

data = fits.open('/gpfs/akrolewski/parity_odd_4pcf/data/random0_DR12v5_CMASS_North.fits.gz')[1].data
data = data[(data['Z'] >= zmin) & (data['Z'] <= zmax) ]
if region == 1:
	data = data[(data['Z'] >= zmin) & (data['Z'] <= zmax) 
		& (data['RA']< (155 - 1./12. * data['DEC']))]
elif region == 2:
	data = data[(data['Z'] >= zmin) & (data['Z'] <= zmax) 
		& (data['RA'] > (165 + 1./12. * data['DEC']))
		& (data['RA'] < (205. - 1./12. * data['DEC'])
		)]
elif region == 3:
	data = data[(data['Z'] >= zmin) & (data['Z'] <= zmax) 
		& (data['RA']> (215. + 1./12. * data['DEC'])
		)]



np.random.seed(42)
data_inds = np.arange(len(data))
np.random.shuffle(data_inds) 

#data = data[:len(data)//2]

#np.random.seed(42)
#data_inds = np.arange(len(data))
#np.random.shuffle(data_inds) 


cosmo = FlatLambdaCDM(H0 = 67.6, Om0 = 0.31, m_nu = [0, 0, 0.06] *u.eV, Ob0 = 0.022/(0.676**2.))

rr = cosmo.comoving_distance(data['Z'][data_inds]).value * 0.676

x = rr * np.sin((90. - data['DEC'][data_inds]) * np.pi/180.) * np.cos(data['RA'][data_inds]*np.pi/180.)
y = rr * np.sin((90. - data['DEC'][data_inds]) * np.pi/180.) * np.sin(data['RA'][data_inds]*np.pi/180.)
z = rr * np.cos((90. - data['DEC'][data_inds]) * np.pi/180.)

weight = data['WEIGHT_FKP'][data_inds]
#if type == 'galaxy':
#elif type == 'random0':
leng = len(data)
for i in range(32):
	#if (i == 0) or (j == 0):
	if i < 31:
		# Note from run_npcf.csh
		# We expect the summed weights to be the same for the data and each random catalog, but the random weights should be negative
		# So we need to re-normalize the weights here and make them negative
		weight_ind = weight[i * leng//32: (i + 1)* leng//32]
	else:
		weight_ind =  weight[i * leng//32: ]
	weight_ind *= sum_wts/np.sum(weight_ind)
	print('sum data weights',sum_wts)
	print('sum random weights',np.sum(weight_ind))
	out = np.array([x[i * leng//32: (i + 1)* leng//32],
			y[i * leng//32: (i + 1)* leng//32],
			z[i * leng//32: (i + 1)* leng//32],
			-1 * weight_ind]).T

	np.savetxt('/gpfs/akrolewski/parity_odd_4pcf/data/random0_DR12v5_CMASS_North_reg%i_cart.ran.%02d' % (region, i),out)

#print('Done with %i' % j, time.time()-t0)
