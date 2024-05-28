import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

hemi = 'South'
type = 'random0'

# Load data and convert into space-separated x, y, z, w in fiducial cosmology

data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data
cosmo = FlatLambdaCDM(H0 = 67.6, Om0 = 0.31, m_nu = [0, 0, 0.06] *u.eV, Ob0 = 0.022/(0.676**2.))

rr = cosmo.comoving_distance(data['Z']).value * 0.676

x = rr * np.sin((90. - data['DEC']) * np.pi/180.) * np.cos(data['RA']*np.pi/180.)
y = rr * np.sin((90. - data['DEC']) * np.pi/180.) * np.sin(data['RA']*np.pi/180.)
z = rr * np.cos((90. - data['DEC']) * np.pi/180.)

if type == 'galaxy':
	weight = data['WEIGHT_SYSTOT'] * (data['WEIGHT_NOZ'] + data['WEIGHT_CP'] - 1.0) * data['WEIGHT_FKP']
elif type == 'random0':
	weight = data['WEIGHT_FKP']

if type == 'galaxy':
	out = np.array([x, y, z, weight]).T
	np.savetxt('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70_cart.data' % (type, hemi),out)
elif type == 'random0':
	leng = len(data)
	for i in range(32):
		if i < 31:
			out = np.array([x[i * leng//32: (i + 1)* leng//32],
				y[i * leng//32: (i + 1)* leng//32],
				z[i * leng//32: (i + 1)* leng//32],
				-1 * weight[i * leng//32: (i + 1)* leng//32]]).T
		else:
			out = np.array([x[i * leng//32: ],
				y[i * leng//32: ],
				z[i * leng//32: ],
				-1 * weight[i * leng//32: ]]).T
		np.savetxt('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70_cart.ran.%02d' % (type, hemi, i),out)
			
