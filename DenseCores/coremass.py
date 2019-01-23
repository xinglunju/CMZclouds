import astropy.constants as const
import astropy.units as u
from astropy.modeling.blackbody import blackbody_nu
import numpy as np

c = const.c.cgs.value
h = const.h.cgs.value
kb = const.k_B.cgs.value
msun = const.M_sun.cgs.value
pc = const.pc.cgs.value
AU = u.AU.cgs.scale
mp = const.m_p.cgs.value
jy = u.Jy.cgs.scale
G = const.G.cgs.value

def coreprop(freq, tdust, flux, dist, kappa, sigmav, tgas, smaj, smin, wmol):
	"""
	Input sigmav is 1-D velocity dispersion, deconvolved with instrument channel width.
	smaj & smin are deconvolved size in arcsec (from casaviewer 2-D fitter)
	"""
	BT = blackbody_nu(freq*u.GHz, tdust*u.K).cgs.value
	d = dist * 1e3 * pc
	# Radius (in unit of pc)
	radius = np.sqrt(smaj*smin)/2/3600/180*np.pi * d / pc
	# Mass
	m = 1e2 * flux * jy * d**2 / BT / kappa / msun
	# Number density
	density = m*msun / (4./3. * np.pi * (radius * pc)**3) / (2.8 * mp)
	# Sigmav of 2.33H
	sigmav = np.sqrt((sigmav*1e5)**2 - kb*tgas/wmol/mp + kb*tgas/2.33/mp) / 1e5
	# Virial
	alpha = 5 * (sigmav * 1e5)**2 * (radius * pc) / G / (m * msun)
	# Virial parameter with B-field, assume 1 mG
	B = 1 * 1e-3 # G
	rho = density * 2.8 * mp # g/cm3
	alfven = B / np.sqrt(4*np.pi*rho)
	alpha_B = 5 * ( (sigmav * 1e5)**2 + (alfven)**2/6.) * (radius * pc) / G / (m * msun)

	return sigmav, radius, m, density, alpha, alfven*1e-5, alpha_B

def sonic(freq, temp, flux, dist, kappa, sigmav, radius, wmol):
	thermalv = np.sqrt(kb * temp / wmol / mp) * 1e-5
	turbv = np.sqrt(sigmav**2 - thermalv**2)
	
	return "%.2f, %.2f, Mach = %.1f" % (thermalv, turbv, turbv/thermalv)

## The upper limit of 'bound mass', using the total 1.3 mm continuum flux in the cloud
print '#######  sigmav radus mass density virial alfvenic virial_B'
print '20kms_allSMA  %.2f %.2f %.2f %.1e %.2f %.2f %.2f' % coreprop(224.9,20,4.5,8.1,0.899,1.14,100,7.22,3.16,29.)
print 'sgrc_allSMA  %.2f %.2f %.2f %.1e %.2f %.2f %.2f' % coreprop(224.9,20,3.47,8.1,0.899,1.14,100,7.22,3.16,29.)

## Masses of cores
print '#######  sigmav radus mass density virial alfvenic virial_B'
print '20kms_c1p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.648*0.5,8.1,0.899,1.14,100,7.22,3.16,29.)
print '20kms_c1p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.192*0.5,8.1,0.899,1.16,100,11.1,2.8,29.)
print '20kms_c1p3  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.113*0.5,8.1,0.899,1.74,100,12.2,3.3,29.)
print '20kms_c2p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.189*0.5,8.1,0.899,2.44,100,10.1,4.1,29.)
print '20kms_c2p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.119*0.5,8.1,0.899,1.32,100,9.1,5.0,29.)
print '20kms_c2p3  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.088*0.5,8.1,0.899,0.98,100,7.0,4.0,29.)
print '20kms_c3p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.208*0.5,8.1,0.899,1.84,100,6.4,3.1,29.)
print '20kms_c3p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.063*0.5,8.1,0.899,1.27,100,9.0,0.5,32.)
print '20kms_c3p3  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.113*0.5,8.1,0.899,1.36,100,9.0,5.1,29.)
print '20kms_c4p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.459,8.1,0.899,1.45,100,7.8,5.4,29.)
print '20kms_c4p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.195,8.1,0.899,1.31,100,14.,3.5,32.)
print '20kms_c4p3  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.104*0.5,8.1,0.899,1.13,100,4.1,1.7,29.)
print '20kms_c4p4  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.076*0.5,8.1,0.899,2.07,100,5.0,2.7,29.)
print '20kms_c4p5  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.035*0.5,8.1,0.899,0.92,100,5.2,3.1,29.)
print '20kms_c4p6  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.038*0.5,8.1,0.899,1.05,100,3.4,2.6,29.)
print '20kms_c5p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.190*0.5,8.1,0.899,1.07,100,7.1,3.9,29.)
print '20kms_c5p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.151*0.5,8.1,0.899,0.89,100,6.6,4.4,29.)


print '#######  sigmav radus mass density virial alfvenic virial_B'
print '50kms_c1p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.143/2,8.1,0.899,5.14,100,4.3,3.1,29.)
print '50kms_c1p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.093/2,8.1,0.899,5.18,100,4.2,2.2,29.)
print '50kms_c2p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.000/2,8.1,0.899,2.57,100,4.0,2.5,29.)
print '50kms_c2p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.000/2,8.1,0.899,1.71,100,2.7,1.8,29.)

print '#######  sigmav radus mass density virial alfvenic virial_B'
print 'G0253_c1p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.061/2,8.1,0.899,2.28,100,5.4,3.7,29.)
print 'G0253_c1p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.114/2,8.1,0.899,2.28,100,7.0,5.2,29.)
print 'G0253_c1p3  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.033/2,8.1,0.899,1.27,100,4.7,3.7,29.)
print 'G0253_c2p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.076,8.1,0.899,2.55,100,10.4,3.2,29.)
print 'G0253_c2p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.083/2,8.1,0.899,1.10,100,5.3,3.7,29.)
print 'G0253_c2p3  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.043/2,8.1,0.899,1.37,100,4.2,2.6,29.)
print 'G0253_c2p4  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.076/2,8.1,0.899,2.56,100,7.7,2.8,29.)
print 'G0253_c2p5  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.041/2,8.1,0.899,1.77,100,7.1,2.9,29.)
print 'G0253_c2p6  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.039/2,8.1,0.899,0.00,100,7.0,4.4,29.)
print 'G0253_c3p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.093/2,8.1,0.899,3.19,100,4.9,2.5,32.)
print 'G0253_c3p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.102/2,8.1,0.899,2.13,100,9.7,4.2,29.)
print 'G0253_c3p3  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.046/2,8.1,0.899,1.50,100,6.8,2.7,29.)
print 'G0253_c3p4  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.158/2,8.1,0.899,3.80,100,10.8,7.6,32.)
print 'G0253_c4p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.093/2,8.1,0.899,1.87,100,9.2,7.1,29.)
print 'G0253_c4p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.127/2,8.1,0.899,1.90,100,17.5,2.4,29.)
print 'G0253_c4p3  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.030/2,8.1,0.899,2.36,100,4.2,1.5,29.)

print '#######  sigmav radus mass density virial alfvenic virial_B'
print 'SgrB1_c1p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.142/2,8.1,0.899,3.14,100,9.2,5.5,32.)
print 'SgrB1_c2p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.221,8.1,0.899,0.61,100,8.4,3.7,29.)
print 'SgrB1_c2p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.289/2,8.1,0.899,1.75,100,6.1,3.6,29.)
print 'SgrB1_c2p3  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.165/2,8.1,0.899,1.71,100,8.3,3.4,29.)
print 'SgrB1_c2p4  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.039,8.1,0.899,2.40,100,8.6,3.6,29.)
print 'SgrB1_c2p5  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.099/2,8.1,0.899,0.00,100,4.6,3.5,29.)
print 'SgrB1_c2p6  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.041,8.1,0.899,1.37,100,6.7,6.4,29.)

print '#######  sigmav radus mass density virial alfvenic virial_B'
print 'SgrC_c1p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.107/2,8.1,0.899,1.08,100,7.3,0.4,32.)
print 'SgrC_c2p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.159/2,8.1,0.899,0.78,100,9.5,6.9,29.)
print 'SgrC_c3p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.205,8.1,0.899,1.50,100,7.1,2.8,29.)
print 'SgrC_c3p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.245/2,8.1,0.899,1.76,100,6.3,2.5,32.)
print 'SgrC_c3p3  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.237/2,8.1,0.899,1.27,100,7.0,1.4,29.)
print 'SgrC_c4p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.457,8.1,0.899,1.71,100,4.2,2.8,29.)
print 'SgrC_c4p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.330,8.1,0.899,1.60,100,4.5,2.9,29.)
print 'SgrC_c5p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.084/2,8.1,0.899,2.11,100,4.1,2.2,29.)
print 'SgrC_c5p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,50,0.059/2,8.1,0.899,2.24,100,3.7,1.7,29.)

print '#######  sigmav radus mass density virial alfvenic virial_B'
print 'SgrD_c1p1  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,25,0.152,2.36,0.899,0.68,100,13.1,8.2,29.)
print 'SgrD_c1p2  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,25,0.144,2.36,0.899,1.59,100,7.9,4.5,29.)
print 'SgrD_c1p3  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,25,0.432/2,2.36,0.899,2.47,100,7.2,2.8,29.)
print 'SgrD_c1p4  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,25,0.095,2.36,0.899,1.75,100,6.8,3.1,29.)
print 'SgrD_c1p5  %.2f %.2f %.2f %.2e %.2f %.2f %.2f' % coreprop(224.9,25,0.157/2,2.36,0.899,1.31,100,6.9,4.0,32.)

