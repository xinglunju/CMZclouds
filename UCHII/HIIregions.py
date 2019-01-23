import numpy as np

import astropy.units as u
import astropy.constants as const

mp=const.m_p.cgs.value
msol = const.M_sun.cgs.value
pc = const.pc.cgs.value
mp=const.m_p.cgs.value
c = const.pc.cgs.value
e = const.e.esu.value
me = const.m_e.cgs.value
kb=const.k_B.cgs.value
h = const.h.cgs.value

### Table from Davies et al. 2011
mass = np.array([6,9,12,15,20,25,30,35,40,50,60,70,80,90,100,110,120])
Qlyc = np.array([43.33,44.76,46.02,47.03,48.00,48.46,48.69,48.90,49.09,49.31,49.43,49.51,49.58,49.65,49.70,49.76,49.81])
###

# Convert cm flux to photon rate:
flux = np.array([0.1483,0.0015,0.4819,0.1291,0.1657,0.1110,0.0032,5e-4,2.3e-3,0.7e-3,5e-3,14.9e-3]) # Jy
T = np.full((len(flux)),1e4) # K

# Equation from Megzer 1974
Q2 = 4.761e48 / 1. * (23.42)**0.1 * (T)**(-0.45) * flux * (8.10)**2

#print np.log10(Q2)
print np.interp(np.log10(Q2), Qlyc, mass)

####################################
print "Sgr D, distance is 2.36 kpc"
flux = [1.0258] # Jy
T = np.full((len(flux)),1e4) # K
Q2 = 4.761e48 / 1. * (23.42)**0.1 * (T)**(-0.45) * flux * (2.36)**2

#print np.log10(Q2)
print np.interp(np.log10(Q2), Qlyc, mass)

