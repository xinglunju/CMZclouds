import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

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
rcParams['xtick.minor.visible'] = 'False'
rcParams['ytick.minor.visible'] = 'True'

# Generate initial figure, scatter plot, and histogram quadrants
fig = plt.figure(111, figsize=(6, 5))

ax = fig.add_subplot(111, position=[0.15,0.12,0.83,0.86])
ax.set_xlabel(r'Number of High-mass Protostars')
ax.set_ylabel(r'$M_\mathrm{cluster}$ ($M_\odot$)')
ax.set_xlim(0, 11)
ax.set_ylim(0., 1300)
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
ax.text(0.05,0.92,r'(b)',color='k',fontsize='x-large',weight='black',transform=ax.transAxes)

N = np.arange(1,11,1)
#M       = np.array([81.82,178.24,274.42,369.70,464.99,561.38,657.35,754.10,849.08,944.48])
#M_sigma = np.array([80.64,112.31,135.13,156.64,173.10,190.11,206.06,217.63,234.74,242.41])
M       = np.array([91.17,198.00,301.32,402.42,502.43,603.33,703.73,803.34,903.00,1000.65])
M_sigma = np.array([81.73,113.75,137.50,158.74,176.34,192.65,209.00,222.93,236.55,247.00])
# Add error bars
sc = ax.errorbar(N, M, yerr=M_sigma, fmt='ko', zorder=10, ecolor='k',elinewidth=0.5,capsize=2.0)

# Linear fit
#m, b = np.polyfit(N, M, 1, w=1/M_sigma)
#ax.plot(N, b+m*N, '--',zorder=10,color='red')
#print m, b

ax.plot(N, 95.8*N, '--', color='blue')

#from scipy.optimize import curve_fit
#def simpleline(N,m):
#	return m*N

#popt,pcov = curve_fit(simpleline,N,M,p0=[96],sigma=M_sigma,absolute_sigma=True)
#ax.plot(N, popt[0]*N, '-',zorder=10,color='red')
#slope = popt[0]
#err   = np.sqrt(np.diag(pcov))[0]

#ax.text(0.35,0.08,r'$M_\mathrm{cluster}=(%.1f\pm%.1f)\times N~M_{\odot}$'%(slope,err),color='red',fontsize='large',weight='black',transform=ax.transAxes)

plt.draw()

# Save figure
fig.savefig('Mcluster_vs_N.pdf')
