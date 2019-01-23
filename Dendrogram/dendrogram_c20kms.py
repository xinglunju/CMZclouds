
## The source code of aplpy is changed to plot this figure
## see line #152 of /Users/xlu/Program/anaconda2/lib/python2.7/site-packages/aplpy/colorbar.py
# : self._colorbar_axes.tick_params(pad=-24,length=7)

import numpy as np
from astrodendro import Dendrogram, pp_catalog
from astropy.io import fits
from astropy import units as u
import aplpy

import matplotlib.pyplot as plt
import matplotlib.colors as colors

from matplotlib import rc, rcParams
from matplotlib.font_manager import fontManager, FontProperties
#rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
rcParams['text.latex.preamble'] = [r'\usepackage{siunitx}',r'\sisetup{detect-all}',r'\usepackage{helvet}',r'\usepackage{sansmath}',r'\sansmath']
rcParams.update({'font.size': 20})
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

####
hdu = fits.open('c20kms.cont.image.fits')[0]

# min_value/min_delta = 4/1 rms
# min_npix = pixels in one beam, pi * bmaj * bmin / 4 / (pixel)^2
rms = 3e-3
beam = 16.63
d = Dendrogram.compute(hdu.data, min_value=4*rms, min_delta=rms, min_npix=beam)

leaves_repro = []
for leaf in d.leaves:
	indices = leaf.indices()
	N_pix = (hdu.data[indices]>5*rms).sum()
	if N_pix > beam:
		leaves_repro.append(leaf)

print 'Number of leaves:', np.size(leaves_repro)

metadata = {}
metadata['data_unit'] = u.Jy / u.beam
metadata['spatial_scale'] =  0.8 * u.arcsec
metadata['beam_major'] =  4.87 * u.arcsec
metadata['beam_minor'] =  2.78 * u.arcsec
cat = pp_catalog(leaves_repro, metadata)
#cat.pprint(show_unit=True, max_lines=200)
fluxes = cat['flux'].data


# Create empty mask. For each leaf we do an 'or' operation with the mask so
# that any pixel corresponding to a leaf is set to True.
mask = np.zeros(hdu.data.shape, dtype=bool)
for leaf in leaves_repro:
	mask = mask | leaf.get_mask()

# Now we create a FITS HDU object to contain this, with the correct header
mask_hdu = fits.PrimaryHDU(mask.astype('short'), hdu.header)

# Load the PB corrected image
#pbim = fits.open('../field3_cont_test.pbcor.image.fits')[0]

# We then use APLpy to make the final plot
fig = plt.figure(figsize=(5.8, 8))

mc = [266.4055, -29.085]
dist = 8.1 # kpc
f1 = aplpy.FITSFigure(hdu,figure=fig,subplot=[0.18,0.08,0.76,0.88])
f1.show_colorscale(cmap='gist_heat_r', vmax=0.196, vmin=0.010, aspect='equal')
f1.show_contour(mask_hdu, colors='blue', linewidths=1.0)
f1.recenter(mc[0], mc[1], width=0.046, height=0.070)
f1.show_contour('../../flux.pb.huge.fits',levels=[0.5],colors='gold',linewidths=2.0,linestyles='dashed',zorder=30)
f1.tick_labels.set_xformat('hh:mm:ss')
f1.tick_labels.set_yformat('dd:mm')
f1.ticks.set_xspacing(4./3600.*15.)
#f1.ticks.set_yspacing(30./3600.)
f1.axis_labels.set_ypad(-27)
f1.set_nan_color('white')
f1.ticks.set_color('black')
f1.ticks.set_minor_frequency(4)

# Scale bar
f1.add_scalebar(0.2/dist/1e3/np.pi*180,label=r'0.2 pc',corner='bottom right',linewidth=2)
f1.scalebar.set_font_size('large')
# Beam
f1.add_beam()
f1.beam.set_frame(True)
f1.beam.set_color('gold')
f1.ticks.set_length(6)
# Color bar
f1.add_colorbar()
f1.colorbar.set_location('top')
f1.colorbar.set_pad(0)
f1.colorbar.set_width(0.1)
f1.colorbar.set_ticks([0.05,0.1,0.15])
#f1.colorbar.set_axis_label_pad(15)
#f1.colorbar.set_axis_label_text('$I_\mathrm{1.3\,mm}$ (mJy/beam)')

f1.add_label(0.05, 0.95, r'20 km\,s$^{-1}$ Cloud', relative=True, weight='black', size='x-large',ha='left')

fig.savefig('c20kms_dendrogram_plot.pdf')
