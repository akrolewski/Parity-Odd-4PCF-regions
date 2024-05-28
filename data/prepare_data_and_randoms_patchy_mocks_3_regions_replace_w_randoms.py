import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import time
import sys

region = int(sys.argv[1])
#region = 2
frac_to_replace_w_randoms = 0.1
rank = int(sys.argv[2])
#rank = 0

t0 = time.time()

hemi = 'NGC'


# Load data and convert into space-separated x, y, z, w in fiducial cosmology

#data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data
zmin = 0.43
zmax = 0.70

for j in range(0,99):
	if rank == j:
		data = np.loadtxt('/gpfs/akrolewski/parity_odd_4pcf/data/patchy_mocks/Patchy-Mocks-DR12%s-COMPSAM_V6C_%04d.dat' % (hemi, j+1))
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

		weight_data = data[:,6] * data[:,7]/(1 + 10000*data[:,4])
		cosmo = FlatLambdaCDM(H0 = 67.6, Om0 = 0.31, m_nu = [0, 0, 0.06] *u.eV, Ob0 = 0.022/(0.676**2.))

		rr = cosmo.comoving_distance(data[:,2]).value * 0.676

		x_data = rr * np.sin((90. - data[:,1]) * np.pi/180.) * np.cos(data[:,0]*np.pi/180.)
		y_data = rr * np.sin((90. - data[:,1]) * np.pi/180.) * np.sin(data[:,0]*np.pi/180.)
		z_data = rr * np.cos((90. - data[:,1]) * np.pi/180.)

		#if type == 'galaxy':
		#weight = data['WEIGHT_SYSTOT'] * (data['WEIGHT_NOZ'] + data['WEIGHT_CP'] - 1.0) * data['WEIGHT_FKP']
		sum_wts = np.sum(weight_data)
		lendata = len(data)

		#out = np.array([x, y, z, weight]).T
		#np.savetxt('/gpfs/akrolewski/parity_odd_4pcf/data/patchy_mocks/Patchy-Mocks-DR12%s_reg%i-COMPSAM_V6C_%04d_cart.data' % (hemi, region, j+1),out)
		#print(j)

		data = np.loadtxt('/gpfs/akrolewski/parity_odd_4pcf/data/patchy_mocks/Patchy-Mocks-Randoms-DR12NGC-COMPSAM_V6C_x50.dat')
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
		np.random.shuffle(data_inds) # Try this for NGC. Didn't do this for SGC
		print('shuffled inds')

		cosmo = FlatLambdaCDM(H0 = 67.6, Om0 = 0.31, m_nu = [0, 0, 0.06] *u.eV, Ob0 = 0.022/(0.676**2.))

		rr = cosmo.comoving_distance(data[:,2][data_inds]).value * 0.676
		print('distance')

		x = rr * np.sin((90. - data[:,1][data_inds]) * np.pi/180.) * np.cos(data[:,0][data_inds]*np.pi/180.)
		y = rr * np.sin((90. - data[:,1][data_inds]) * np.pi/180.) * np.sin(data[:,0][data_inds]*np.pi/180.)
		z = rr * np.cos((90. - data[:,1][data_inds]) * np.pi/180.)

		weight = data[:,5][data_inds] * data[:,6][data_inds]/(1 + 10000*data[:,3][data_inds])

		x_data = np.concatenate((x_data, x[:int(frac_to_replace_w_randoms * lendata)]))
		y_data = np.concatenate((y_data, y[:int(frac_to_replace_w_randoms * lendata)]))
		z_data = np.concatenate((z_data, z[:int(frac_to_replace_w_randoms * lendata)]))

		weight_data = np.concatenate((weight_data, weight[:int(frac_to_replace_w_randoms * lendata)]))

		#ra_data = np.concatenate((ra_data, data[:,0][:int(frac_to_replace_w_randoms * lendata)]))
		#dec_data = np.concatenate((dec_data, data[:,1][:int(frac_to_replace_w_randoms * lendata)]))
		#redshift_data  = np.concatenate((redshift_data, data[:,2][:int(frac_to_replace_w_randoms * lendata)]))
		#weight_data = np.concatenate((weight_data, np.ones_like(data['RA'][:int(frac_to_replace_w_randoms * lendata)])))
		#weight_fkp_data = np.concatenate((weight_fkp_data, data['WEIGHT_FKP'][:int(frac_to_replace_w_randoms * lendata)]))
		#nz_data  = np.concatenate((nz_data, data['NZ'][:int(frac_to_replace_w_randoms * lendata)]))
		print('made arrays')


		sum_wts = np.sum(weight_data)


		out = np.array([x_data, y_data, z_data, weight_data]).T
		np.savetxt('/gpfs/akrolewski/parity_odd_4pcf/data/patchy_mocks/Patchy-Mocks-DR12%s_reg%i_replace_%.2f_w_random-COMPSAM_V6C_%04d_cart.data' % (hemi, region, frac_to_replace_w_randoms, j+1),out)

		# from astropy.table import Table
		# tdat = Table(np.array([ra_data.astype('>f8'), dec_data.astype('>f8'), redshift_data.astype('>f8'), weight_data.astype('>f8'), weight_fkp_data.astype('>f8'), nz_data.astype('>f8')]).T, 
		# 	names=('RA','DEC','Z','WEIGHT','WEIGHT_FKP','NZ'), dtype=('>f8','>f8','>f8','>f8','>f8','>f8'))
		# #t.write('galaxy_DR12v5_CMASSLOWZTOT_%s_reg%i_replace_%.2f_w_random_zcut_0.43_0.70.fits' % (hemi, region, frac_to_replace_w_randoms), overwrite=True)

		x = x[int(frac_to_replace_w_randoms * lendata):]
		y = y[int(frac_to_replace_w_randoms * lendata):]
		z = z[int(frac_to_replace_w_randoms * lendata):]

		weight = weight[int(frac_to_replace_w_randoms * lendata):]

		# tran = Table(np.array([data['RA'][data_inds][int(frac_to_replace_w_randoms * lendata):],
		# 	data['DEC'][data_inds][int(frac_to_replace_w_randoms * lendata):],
		# 	data['Z'][data_inds][int(frac_to_replace_w_randoms * lendata):],
		# 	np.ones_like(data['Z'][data_inds][int(frac_to_replace_w_randoms * lendata):]),
		# 	data['WEIGHT_FKP'][data_inds][int(frac_to_replace_w_randoms * lendata):],
		# 	data['NZ'][data_inds][int(frac_to_replace_w_randoms * lendata):]]).T, names=('RA','DEC','Z','WEIGHT','WEIGHT_FKP','NZ'))
		# #t.write('random0_DR12v5_CMASSLOWZTOT_%s_reg%i_replace_%.2f_w_random_zcut_0.43_0.70.fits' % (hemi, region, frac_to_replace_w_randoms))


		#if ttype == 'galaxy':
		#elif ttype == 'random0':
		leng = len(weight)
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

				np.savetxt('/gpfs/akrolewski/parity_odd_4pcf/data/patchy_mocks/Patchy-Mocks-DR12%s_reg%i_replace_%.2f_w_random-COMPSAM_V6C_x50_%04d_cart.ran.%02d' % (hemi, region, frac_to_replace_w_randoms, j+1, i),out)
