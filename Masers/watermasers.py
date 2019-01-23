import numpy as np
import astropy.units as u
import astropy.constants as const

# Maser fluxes
I_table_cmz = [1090,491,946,4320,852,7153,418,242,49,100,337,3738,411,153,327,2762,94,360,95,349,1432,263,2790,17,36,541,372,73,188,21,945,364,24,469,587,40,301,43,96,9037,2495,105,373,2760,42370,12700,596,8152,8360]
# in mJy * km/s
Imaser = np.array(I_table_cmz) / 1e3 # in Jy
bdwidth = 1.0/3e5*22.235*1e9 # in Hz
dist = 8100 * const.pc.cgs.value # in cm
Lmaser_cmz = Imaser * bdwidth * 4 * np.pi * (dist)**2 * 1e-23 / const.L_sun.cgs.value # in Lsol

# Masers in Sgr D
I_table_sgrd = [11450,110,5144,40,86,7302,576]
Imaser = np.array(I_table_sgrd) / 1e3 # in Jy
dist = 2360 * const.pc.cgs.value # in cm
Lmaser_sgrd = Imaser * bdwidth * 4 * np.pi * (dist)**2 * 1e-23 / const.L_sun.cgs.value # in Lsol

Lmaser = np.concatenate([Lmaser_cmz, Lmaser_sgrd])

#print "Maser luminosities (10^-7 Lsol)", Lmaser*1e7

# Correlation from Urquhart 2011
Lstar = (Lmaser / 7.1e-12)**(1/1.47) # in Lsol
print np.c_[np.concatenate([I_table_cmz,I_table_sgrd]), Lmaser*1e7, np.log10(Lstar)]

# The M-L relation from Davies 2011
mass = np.array([6,9,12,15,20,25,30,35,40,50,60,70,80,90,100,110,120])
lum =  np.array([3.01,3.57,3.96,4.26,4.61,4.86,5.02,5.18,5.34,5.52,5.70,5.81,5.92,6.02,6.09,6.16,6.23])

## 20kms W1, W2, W3, W5, W15
l = [3.416, 3.610, 4.059, 4.208, 3.926]
print np.interp(l, lum, mass)

# G0.253 W2
l = [3.445]
print np.interp(l, lum, mass)

# SgrC W2, W7, W11, W14
l = [3.469, 4.277, 3.926, 3.473]
print np.interp(l, lum, mass)
