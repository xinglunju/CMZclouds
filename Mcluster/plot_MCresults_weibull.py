import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
#import scipy.special as sp

from matplotlib import rc, rcParams
from matplotlib.font_manager import fontManager, FontProperties
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
rcParams['text.latex.preamble'] = [ r'\usepackage{siunitx}', r'\sisetup{detect-all}', r'\usepackage{helvet}', r'\usepackage{sansmath}', r'\sansmath']
rcParams.update({'font.size': 14})
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['xtick.top'] = 'True'
rcParams['ytick.right'] = 'True'
rcParams['xtick.minor.visible'] = 'True'
rcParams['ytick.minor.visible'] = 'True'

m, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10 = np.loadtxt('HMdet_10_dr_1.0_log.dat', usecols = (0,1,2,3,4,5,6,7,8,9,10), unpack = True)

# Generate initial figure, scatter plot, and histogram quadrants
fig = plt.figure(111, figsize=(6, 5))

ax = fig.add_subplot(111, position=[0.15,0.12,0.83,0.86])
ax.set_xlabel(r'$M_\mathrm{cluster}$ ($M_\odot$)')
ax.set_ylabel(r'Normalized Probability')
ax.set_xlim(8, 3000)
ax.set_ylim(0., 1.1)
ax.set_xscale('log')
#ax.set_yscale('log')

l1 /= l1.max()
l2 /= l2.max()
l3 /= l3.max()
l4 /= l4.max()
l5 /= l5.max()
l6 /= l6.max()
l7 /= l7.max()
l8 /= l8.max()
l9 /= l9.max()
l10 /= l10.max()
ax.plot(m, l1, 'bo')
#ax.plot(m, l2, 'ro')
#ax.plot(m, l3, 'go')
ax.plot(m, l4, 'ko')
#ax.plot(m, l5, 'cx')
#ax.plot(m, l6, 'm+')
ax.plot(m, l7, '*', color='gray')
#ax.plot(m, l8, 's', color='orange')
#ax.plot(m, l9, '^', color='b')
ax.plot(m, l10, '^', color='r')

ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.tick_params(axis='both',which='both',zorder=20)
ax.grid(color='gray', linestyle='-', linewidth=0.5, axis='x', which='both')
ax.text(0.05,0.92,r'(a)',color='k',fontsize='x-large',weight='black',transform=ax.transAxes)

#ax.axvline(x=30,linestyle='--',color='black',linewidth=0.5)
#ax.axhline(y=2,linestyle='--',color='black',linewidth=0.5)
#ax.axhline(y=1,linestyle='--',color='black',linewidth=0.5)

#def weibull(x,x0,sigma):
#	return np.exp(-0.5*(x-x0)**2/sigma**2)

## Weibull function
def weibull(x,lam,k):
	wbfunc = k/lam * (x/lam)**(k-1) * np.exp(-(x/lam)**k)
	return wbfunc/wbfunc.max()

#popt,pcov = curve_fit(weibull,np.log10(m),l1,p0=[1.9,1.8])
#print(popt)
#print 'Peak is at:', 10**(popt[0])
#fwhm = 10**(popt[0]+np.sqrt(2*np.log(2))*np.abs(popt[1])) - 10**(popt[0]-np.sqrt(2*np.log(2))*np.abs(popt[1]))
#print 'FWHM is:', fwhm
#print 'Sigma is:', fwhm / 2.355
#plt.plot(m,weibull(np.log10(m),*popt),'b-')

popt,pcov = curve_fit(weibull,np.log10(m),l1,p0=[2.0,5.0])
print(popt)
print 'Peak is at:', 10**(popt[0]*(1.-1./popt[1])**(1./popt[1]))
l1 /= l1.sum()
mean = (m*l1).sum()
print 'RMS is:', np.sqrt(((m - mean)**2*l1).sum())
plt.plot(m,weibull(np.log10(m),*popt),'b-')

popt,pcov = curve_fit(weibull,np.log10(m),l2,p0=[2.0,5])
print(popt)
print 'Peak is at:', 10**(popt[0]*(1.-1./popt[1])**(1./popt[1]))
l2 /= l2.sum()
mean = (m*l2).sum()
print 'RMS is:', np.sqrt(((m - mean)**2*l2).sum())
plt.plot(m,weibull(np.log10(m),*popt),'r-')

popt,pcov = curve_fit(weibull,np.log10(m),l3,p0=[2.1,5])
print(popt)
print 'Peak is at:', 10**(popt[0]*(1.-1./popt[1])**(1./popt[1]))
l3 /= l3.sum()
mean = (m*l3).sum()
print 'RMS is:', np.sqrt(((m - mean)**2*l3).sum())
plt.plot(m,weibull(np.log10(m),*popt),'g-')

popt,pcov = curve_fit(weibull,np.log10(m),l4,p0=[2.2,5])
print(popt)
print 'Peak is at:', 10**(popt[0]*(1.-1./popt[1])**(1./popt[1]))
l4 /= l4.sum()
mean = (m*l4).sum()
print 'RMS is:', np.sqrt(((m - mean)**2*l4).sum())
plt.plot(m,weibull(np.log10(m),*popt),'k-')

popt,pcov = curve_fit(weibull,np.log10(m),l5,p0=[2.2,5])
print(popt)
print 'Peak is at:', 10**(popt[0]*(1.-1./popt[1])**(1./popt[1]))
l5 /= l5.sum()
mean = (m*l5).sum()
print 'RMS is:', np.sqrt(((m - mean)**2*l5).sum())
plt.plot(m,weibull(np.log10(m),*popt),'c-')

popt,pcov = curve_fit(weibull,np.log10(m),l6,p0=[2.3,5])
print(popt)
print 'Peak is at:', 10**(popt[0]*(1.-1./popt[1])**(1./popt[1]))
l6 /= l6.sum()
mean = (m*l6).sum()
print 'RMS is:', np.sqrt(((m - mean)**2*l6).sum())
plt.plot(m,weibull(np.log10(m),*popt),'m-')

popt,pcov = curve_fit(weibull,np.log10(m),l7,p0=[2.4,5])
print(popt)
print 'Peak is at:', 10**(popt[0]*(1.-1./popt[1])**(1./popt[1]))
l7 /= l7.sum()
mean = (m*l7).sum()
print 'RMS is:', np.sqrt(((m - mean)**2*l7).sum())
plt.plot(m,weibull(np.log10(m),*popt),'-',color='gray')

popt,pcov = curve_fit(weibull,np.log10(m),l8,p0=[2.5,5])
print(popt)
print 'Peak is at:', 10**(popt[0]*(1.-1./popt[1])**(1./popt[1]))
l8 /= l8.sum()
mean = (m*l8).sum()
print 'RMS is:', np.sqrt(((m - mean)**2*l8).sum())
plt.plot(m,weibull(np.log10(m),*popt),'-',color='orange')

popt,pcov = curve_fit(weibull,np.log10(m),l9,p0=[2.6,5])
print(popt)
print 'Peak is at:', 10**(popt[0]*(1.-1./popt[1])**(1./popt[1]))
l9 /= l9.sum()
mean = (m*l9).sum()
print 'RMS is:', np.sqrt(((m - mean)**2*l9).sum())
plt.plot(m,weibull(np.log10(m),*popt),'-',color='b')

popt,pcov = curve_fit(weibull,np.log10(m),l10,p0=[2.7,5])
print(popt)
print 'Peak is at:', 10**(popt[0]*(1.-1./popt[1])**(1./popt[1]))
l10 /= l10.sum()
mean = (m*l10).sum()
print 'RMS is:', np.sqrt(((m - mean)**2*l10).sum())
plt.plot(m,weibull(np.log10(m),*popt),'--',color='r')

plt.draw()
plt.close()

# Save figure
fig.savefig('MCresults.pdf')
