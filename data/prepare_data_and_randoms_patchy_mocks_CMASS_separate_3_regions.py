import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import time
import sys
import os

rank = int(sys.argv[1])
region = int(sys.argv[2])
#rank = 0

t0 = time.time()

hemi = 'N'


# Load data and convert into space-separated x, y, z, w in fiducial cosmology

#data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data
zmin = 0.43
zmax = 0.70

for j in range(500):
	if rank == j%64:
		#print(j)
		#j = 0
		#j = 0
		#if not os.path.exists('/gpfs/akrolewski/parity_odd_4pcf/data/patchy_mocks/Patchy-Mocks-DR12CMASS/cart_files_default/Patchy-Mocks-DR12CMASS-%s-V6C-Portsmouth-mass_%04d_cart.ran.%02d' % (hemi, j+1, 0)):
		data = np.loadtxt('/gpfs/akrolewski/parity_odd_4pcf/data/patchy_mocks/Patchy-Mocks-DR12CMASS/radecz_files/Patchy-Mocks-DR12CMASS-%s-V6C-Portsmouth-mass_%04d.dat' % (hemi, j+1))
		if region == 1:
			data = data[(data[:,2] >= zmin) & (data[:,2] <= zmax) & (data[:,6] * data[:,7] > 0)
				& (data[:,0] < (155 - 1./12. * data[:,1]))]
		elif region == 2:
			data = data[(data[:,2] >= zmin) & (data[:,2] <= zmax) & (data[:,6] * data[:,7] > 0)
				& (data[:,0] > (165 + 1./12. * data[:,1]))
				& (data[:,0] < (205. - 1./12. * data[:,1])
				)]
		elif region == 3:
			data = data[(data[:,2] >= zmin) & (data[:,2] <= zmax) & (data[:,6] * data[:,7] > 0)
				& (data[:,0] > (215. + 1./12. * data[:,1])
				)]

		weight = data[:,6] * data[:,7]/(1 + 10000*data[:,4])
		cosmo = FlatLambdaCDM(H0 = 67.6, Om0 = 0.31, m_nu = [0, 0, 0.06] *u.eV, Ob0 = 0.022/(0.676**2.))

		rr = cosmo.comoving_distance(data[:,2]).value * 0.676

		x = rr * np.sin((90. - data[:,1]) * np.pi/180.) * np.cos(data[:,0]*np.pi/180.)
		y = rr * np.sin((90. - data[:,1]) * np.pi/180.) * np.sin(data[:,0]*np.pi/180.)
		z = rr * np.cos((90. - data[:,1]) * np.pi/180.)

		#if type == 'galaxy':
		#weight = data['WEIGHT_SYSTOT'] * (data['WEIGHT_NOZ'] + data['WEIGHT_CP'] - 1.0) * data['WEIGHT_FKP']
		sum_wts = np.sum(weight)

		out = np.array([x, y, z, weight]).T
		np.savetxt('/gpfs/akrolewski/parity_odd_4pcf/data/patchy_mocks/Patchy-Mocks-DR12CMASS/reg%i/Patchy-Mocks-DR12CMASS-%s-reg%i-V6C-Portsmouth-mass_%04d_cart.data' % (region, hemi, region, j+1),out)
		print(j)

		data = np.loadtxt('/gpfs/akrolewski/parity_odd_4pcf/data/patchy_mocks/Random-DR12CMASS-%s-V6C-x100.dat' % hemi)
		if region == 1:
			data = data[(data[:,2] >= zmin) & (data[:,2] <= zmax) & (data[:,5] * data[:,6] > 0)
				& (data[:,0] < (155 - 1./12. * data[:,1]))]
		elif region == 2:
			data = data[(data[:,2] >= zmin) & (data[:,2] <= zmax) & (data[:,5] * data[:,6] > 0)
				& (data[:,0] > (165 + 1./12. * data[:,1]))
				& (data[:,0] < (205. - 1./12. * data[:,1]))
				]
		elif region == 3:
			data = data[(data[:,2] >= zmin) & (data[:,2] <= zmax) & (data[:,5] * data[:,6] > 0)
				& (data[:,0] > (215. + 1./12. * data[:,1]))
				]

		np.random.seed(42)
		data_inds = np.arange(len(data))
		np.random.shuffle(data_inds) 

		data = data[:len(data)//2]
	
		np.random.seed(42)
		data_inds = np.arange(len(data))
		np.random.shuffle(data_inds) 


		cosmo = FlatLambdaCDM(H0 = 67.6, Om0 = 0.31, m_nu = [0, 0, 0.06] *u.eV, Ob0 = 0.022/(0.676**2.))

		rr = cosmo.comoving_distance(data[:,2][data_inds]).value * 0.676

		x = rr * np.sin((90. - data[:,1][data_inds]) * np.pi/180.) * np.cos(data[:,0][data_inds]*np.pi/180.)
		y = rr * np.sin((90. - data[:,1][data_inds]) * np.pi/180.) * np.sin(data[:,0][data_inds]*np.pi/180.)
		z = rr * np.cos((90. - data[:,1][data_inds]) * np.pi/180.)

		weight = data[:,5][data_inds] * data[:,6][data_inds]/(1 + 10000*data[:,3][data_inds])

		#if type == 'galaxy':
		#elif type == 'random0':
		leng = len(data)
		for i in range(32):
			if (i == 0) or (j == 0):
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

				np.savetxt('/gpfs/akrolewski/parity_odd_4pcf/data/patchy_mocks/Patchy-Mocks-DR12CMASS/reg%i/Patchy-Mocks-DR12CMASS-%s-reg%i-V6C-Portsmouth-mass_%04d_cart.ran.%02d' % (region, hemi, region, j+1, i),out)

		print('Done with %i' % j, time.time()-t0)
