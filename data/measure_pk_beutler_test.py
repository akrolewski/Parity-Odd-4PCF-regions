import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import healpy as hp
import matplotlib.pyplot as plt
from nbodykit.source.catalog import ArrayCatalog
from nbodykit.lab import *
import sys

#from healpy.newvisufunc import projview, newprojplot
#from mpi4py import MPI

#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#rank = 1
rank = int(sys.argv[1])

zmin = 0.5
zmax = 0.75


hemi = 'North'
type = 'galaxy'

# the fiducial BOSS DR12 cosmology
cosmo = cosmology.Cosmology(h=0.676).match(Omega0_m=0.31)

cosmo_patchy = cosmology.Cosmology(h=0.6777,Omega0_b=0.048, m_ncdm = []).match(Omega0_m=0.307115).match(sigma8=0.8288)

#nside = 128

# Load data and convert into space-separated x, y, z, w in fiducial cosmology
#map_dat = np.zeros(12*nside**2)
#for i in range(2048):
#	if i%32 == rank:
i = rank + 2
#data = fits.open('%s_DR12v5_CMASSLOWZTOT_%s_zcut_0.43_0.70.fits' % (type, hemi))[1].data
#data = data[(data['z'] >= zmin) & (data['z'] <= zmax)]
#data_wt = data['weight_fkp'] * data['weight_systot'] * (data['weight_cp'] + data['weight_noz']- 1)
data = np.loadtxt('patchy_mocks/Patchy-Mocks-DR12NGC-COMPSAM_V6C_%04d.dat' % (i + 1))
data = data[(data[:,2] >= zmin) & (data[:,2] <= zmax) & (data[:,6] * data[:,7] != 0)]
data_wt = data[:,6] * data[:,7]

nbar_rescaling = cosmo_patchy.comoving_distance(data[:,2])**2./cosmo.comoving_distance(data[:,2])**2. * cosmo.efunc(data[:,2])/cosmo_patchy.efunc(data[:,2])
nbar = data[:,4]# * nbar_rescaling

data_fkp_wt = 1/(1 + 10000*nbar)

data_cat = ArrayCatalog({'RA' : data[:,0], 'DEC' : data[:,1], 'Z': data[:,2], 'Weight': data_wt, 'Weight_FKP': data_fkp_wt, 'NZ': nbar})
		
random = np.loadtxt('patchy_mocks/Patchy-Mocks-Randoms-DR12NGC-COMPSAM_V6C_x50.dat')
random = random[(random[:,2] >= zmin) & (random[:,2] <= zmax) & (random[:,5] * random[:,6] != 0)]
random_wt = random[:,5] * random[:,6]
nbar_rescaling = cosmo_patchy.comoving_distance(random[:,2])**2./cosmo.comoving_distance(random[:,2])**2. * cosmo.efunc(random[:,2])/cosmo_patchy.efunc(random[:,2])
nbar = random[:,3] #* nbar_rescaling

random_fkp_wt = 1/(1 + 10000*nbar)


random_cat = ArrayCatalog({'RA' : random[:,0], 'DEC' : random[:,1], 'Z': random[:,2], 'Weight': random_wt, 'Weight_FKP': random_fkp_wt, 'NZ': nbar })

# add Cartesian position column
data_cat['Position'] = transform.SkyToCartesian(data_cat['RA'], data_cat['DEC'], data_cat['Z'], cosmo=cosmo)
random_cat['Position'] = transform.SkyToCartesian(random_cat['RA'], random_cat['DEC'], random_cat['Z'], cosmo=cosmo)

# combine the data and randoms into a single catalog
fkp = FKPCatalog(data_cat, random_cat)

f = open('patchy_mocks/z3/ps1D_patchyDR12_ngc_combined_bianchi_z3_%i_COMPnbar_TIC_340_650_360_120.dat' % (i+1)).readlines()
Lx = float(f[13].split(' = ')[1].split()[0])
Ly = float(f[13].split(' = ')[2].split()[0])
Lz = float(f[13].split(' = ')[3].split()[0])
nx = float(f[14].split(' = ')[1].split()[0])
ny = float(f[14].split(' = ')[2].split()[0])
nz = float(f[14].split(' = ')[3].split()[0])
xoff = float(f[15].split(' = ')[1].split()[0])
yoff = float(f[15].split(' = ')[2].split()[0])
zoff = float(f[15].split(' = ')[3].split()[0])
mesh = fkp.to_mesh(BoxSize=[Lx, Ly, Lz], BoxCenter=[xoff, yoff, zoff], Nmesh=[nx, ny, nz], nbar='NZ', fkp_weight='Weight_FKP', comp_weight='Weight', window='tsc')


# compute the multipoles
r = ConvolvedFFTPower(mesh, poles=[0,2,4], dk=0.005, kmin=0.)

for key in r.attrs:
    print("%s = %s" % (key, str(r.attrs[key])))
    
# Nbodykit shot noise calculation is wrong
fc = 0.5
alpha_prime = np.sum(random_wt * random_fkp_wt) / np.sum(data_wt * data_fkp_wt)
shot_noise = (np.sum(fc * data_wt * data[:,6] * data_fkp_wt ** 2 + (1 - fc) * data_wt**2 * data_fkp_wt ** 2)
+ np.sum(1./alpha_prime**2 * random_fkp_wt**2))

def shotnoise(r, alpha):
	r"""
	Compute the power spectrum shot noise, using either the
	``data`` or ``randoms`` source.
	This computes:
	.. math::
		S = \sum (w_\mathrm{comp} w_\mathrm{fkp})^2
	References
	----------
	see Eq. 15 of Beutler et al. 2014, "The clustering of galaxies in the
	SDSS-III Baryon Oscillation Spectroscopic Survey: testing gravity with redshift
	space distortions using the power spectrum multipoles"
	"""
	#if 'shotnoise' not in self.attrs:

	Pshot = 0
	for name in ['data', 'randoms']:

		# the selection (same for first/second)
		sel = r.first.source.compute(r.first.source[name][r.first.selection])

		# selected first/second meshes for "name" (data or randoms)
		first = r.first.source[name][sel]
		second = r.second.source[name][sel]

		# completeness weights (assumed same for first/second)
		comp_weight = first[r.first.comp_weight]

		# different weights allowed for first and second mesh
		fkp_weight1 = first[r.first.fkp_weight]
		if r.first is r.second:
			fkp_weight2 = fkp_weight1
		else:
			fkp_weight2 = second[r.second.fkp_weight]

		S = (comp_weight**2*fkp_weight1*fkp_weight2).sum()
		if name == 'randoms':
			S *= alpha**2
		Pshot += S # add to total

	# reduce sum across all ranks
	Pshot = r.comm.allreduce(first.compute(Pshot))

	# divide by normalization from randoms
	return Pshot / r.attrs['randoms.norm']
	
#print(shotnoise(r, 1./alpha_prime))


def shotnoise(r, alpha, weight_collision, weight_veto, weight_fkp, weight_ran_fkp, fc):
	r"""
	Compute the power spectrum shot noise, using either the
	``data`` or ``randoms`` source.
	This computes:
	.. math::
		S = \sum (w_\mathrm{comp} w_\mathrm{fkp})^2
	References
	----------
	see Eq. 15 of Beutler et al. 2014, "The clustering of galaxies in the
	SDSS-III Baryon Oscillation Spectroscopic Survey: testing gravity with redshift
	space distortions using the power spectrum multipoles"
	"""
	#if 'shotnoise' not in self.attrs:

	Pshot = 0
	for name in ['data', 'randoms']:

		# the selection (same for first/second)
		sel = r.first.source.compute(r.first.source[name][r.first.selection])

		# selected first/second meshes for "name" (data or randoms)
		first = r.first.source[name][sel]
		second = r.second.source[name][sel]

		# completeness weights (assumed same for first/second)
		comp_weight = first[r.first.comp_weight]

		# different weights allowed for first and second mesh
		fkp_weight1 = first[r.first.fkp_weight]
		if r.first is r.second:
			fkp_weight2 = fkp_weight1
		else:
			fkp_weight2 = second[r.second.fkp_weight]

		S = (comp_weight**2*fkp_weight1*fkp_weight2).sum()
		if name == 'randoms':
			S *= alpha**2
		Pshot += S # add to total

	# reduce sum across all ranks
	Pshot = r.comm.allreduce(first.compute(Pshot))
	S1 = (np.sum(fc * weight_collision * weight_veto * weight_fkp**2.)
		+ np.sum((1-fc) * weight_collision**2 * weight_fkp**2.))
	S2 = np.sum(alpha**2 * weight_ran_fkp**2)

	# divide by normalization from randoms
	#return Pshot / r.attrs['randoms.norm']
	return (S1 + S2)/r.attrs['randoms.norm']
	
SN = (shotnoise(r, 1./alpha_prime, data[:,7], data[:,6], data_fkp_wt, random_fkp_wt, 0.5))

f = open('patchy_mocks/z3_AK/%i.txt' % (i+1),'w')
f.write('# Shot noise: %.5f\n' % SN)
poles = r.poles
for i in range(len(poles['power_0'].real)):
	f.write('%.5e %.5e %.5e\n' % (poles['power_0'].real[i] - SN, poles['power_2'].real[i], poles['power_4'].real[i]))
f.close()
