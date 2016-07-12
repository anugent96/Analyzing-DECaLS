"""
Makes histograms of rchi2 values for the different types of objects. Plots based on filter used (g, r, z). Also prints out all of the 
rchi2 values for unmasked cosmic rays.
"""

import sys
image = str(sys.argv[1])

from astropy.io import fits
tractor = fits.open(image)
tbl = tractor[1].data

RA = tbl.field('ra') # RA values
DEC = tbl.field('dec') # DEC values
rchi2 = tbl.field('decam_rchi2') # decam rchi2
type1 = tbl.field('type') # type of object (PSF, SIMP, COMP, EXP, DEV)
mask = tbl.field('decam_anymask') # masking values
d_nobs = tbl.field('decam_nobs') # Decam number of exposures
fD = tbl.field('decam_flux') # model flux, decam
fW = tbl.field('wise_flux')
ivd = tbl.field('decam_flux_ivar') # inverse variance (1/sigma squared)
ivw = tbl.field('wise_flux_ivar')

g_rchi2 = rchi2[:,1] # g component of rchi2
r_rchi2 = rchi2[:,2] # r component of rchi2
z_rchi2 = rchi2[:,4] # z component of rchi2

g = fD[:,1] # g flux
r = fD[:,2] # r flux
z = fD[:,4] # z flux
w1 = fW[:,0] # w1 flux

znobs = d_nobs[:,4] # number of observations in z

aMask = [sum(x) for x in mask] # sum of all masking values for 1 object

import matplotlib.pyplot as plt
import math
import itertools
import numpy as np

g_sig2 = np.array(ivd[:,1])
r_sig2 = np.array(ivd[:,2])
z_sig2 = np.array(ivd[:,4])
w1_sig2 = np.array(ivw[:,0])
z_rchi2 = rchi2[:,4]

"""
Flux times the squareroot of its corresponding inverse variance. Gives a better value for flux that can be used to find fluxes that
are significant and not significant (ie: objects that exist and don't exist)
"""
g_test = g*(g_sig2**0.5) 
r_test = r*(r_sig2**0.5)
z_test = z*(z_sig2**0.5)
w1_test = w1*(w1_sig2**0.5)

"""
rchi2 values for possible CRs. Looking for when z-flux is significantly positive and all other fluxes
(g, r, w1) are insignificant.
"""
g1_rchi2 = [x for (x, mask1, mask2) in zip (g_rchi2, w1_test, z_test) if mask1 < 1 and mask2 > 2]
r1_rchi2 = [x for (x, mask1, mask2) in zip (r_rchi2, w1_test, z_test) if mask1 < 1 and mask2 > 2]
z1_rchi2 = [x for (x, mask1, mask2) in zip (z_rchi2, w1_test, z_test) if mask1 < 1 and mask2 > 2]
g1 = [x for (x, mask1, mask2) in zip(g_test, w1_test, z_test) if mask1 < 1 and mask2 > 2]
r1 = [x for (x, mask1, mask2) in zip(r_test, w1_test, z_test) if mask1 < 1 and mask2 > 2]

g2_rchi2 = [x for (x, mask1, mask2) in zip (g1_rchi2, g1, r1) if mask1 < 1 and mask2 < 1]
r2_rchi2 = [x for (x, mask1, mask2) in zip (r1_rchi2, g1, r1) if mask1 < 1 and mask2 < 1]
z2_rchi2 = [x for (x, mask1, mask2) in zip (z1_rchi2, g1, r1) if mask1 < 1 and mask2 < 1]


# rchi2 values in each filter for PSF 
g_psf = [x for (x,y,w) in zip(g_rchi2, type1, g_rchi2) if y == 'PSF' and not math.isnan(w)] 
r_psf = [x for (x,y,w) in zip(r_rchi2, type1, r_rchi2) if y == 'PSF' and not math.isnan(w)]
z_psf = [x for (x,y,w) in zip(z_rchi2, type1, z_rchi2) if y == 'PSF' and not math.isnan(w)] 

# rchi2 values in each filter for possible unmasked CRs
g2_rchi2_psf = [x for (x,y,w) in zip(g2_rchi2, aMask, g2_rchi2) if y == 0 and not math.isnan(w)]
r2_rchi2_psf = [x for (x,y,w) in zip(r2_rchi2, aMask, r2_rchi2) if y == 0 and not math.isnan(w)]
z2_rchi2_psf = [x for (x,y,w) in zip(z2_rchi2, aMask, z2_rchi2) if y == 0 and not math.isnan(w)]

# rchi2 values in each filter for EXP
g_exp = [x for (x,y,w) in zip(g_rchi2, type1, g_rchi2) if y == 'EXP' and not math.isnan(w)] 
r_exp = [x for (x,y,w) in zip(r_rchi2, type1, r_rchi2) if y == 'EXP' and not math.isnan(w)]
z_exp = [x for (x,y,w) in zip(z_rchi2, type1, z_rchi2) if y == 'EXP' and not math.isnan(w)] 

# rchi2 values in each filter for DEV
g_dev = [x for (x,y,w) in zip(g_rchi2, type1, g_rchi2) if y == 'DEV' and not math.isnan(w)] 
r_dev = [x for (x,y,w) in zip(r_rchi2, type1, r_rchi2) if y == 'DEV' and not math.isnan(w)]
z_dev = [x for (x,y,w) in zip(z_rchi2, type1, z_rchi2) if y == 'DEV' and not math.isnan(w)] 

# rchi2 values in each filter for COMP 
g_comp = [x for (x,y,w) in zip(g_rchi2, type1, g_rchi2) if y == 'COMP' and not math.isnan(w)] 
r_comp = [x for (x,y,w) in zip(r_rchi2, type1, r_rchi2) if y == 'COMP' and not math.isnan(w)]
z_comp = [x for (x,y,w) in zip(z_rchi2, type1, z_rchi2) if y == 'COMP' and not math.isnan(w)] 

# g_psf histogram
plt.hist(g_psf, color='k', range=(0,100), bins=150, label='PSF', histtype='stepfilled')
plt.hist(g_exp, color='c', range=(0,100), bins=150, alpha=0.7, label='EXP', histtype='stepfilled')
plt.hist(g_dev, color='y', range=(0,100), bins=150, alpha=0.7, label='DEV', histtype='stepfilled')
plt.hist(g_comp, color='m', range=(0,100), bins=150, alpha=0.7, label= 'COMP', histtype='stepfilled')
plt.hist(g2_rchi2_psf, color='r', range=(0,100), bins=150, alpha=0.9, label='Unmasked CRs', histtype='stepfilled')
plt.title('RChi2 Values (g-filter)')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()

# r_psf histogram
plt.hist(r_psf, color='k', range=(0,100), bins=150, label='PSF', histtype='stepfilled')
plt.hist(r_exp, color='c', range=(0,100), bins=150, alpha=0.7, label='EXP', histtype='stepfilled')
plt.hist(r_dev, color='y', range=(0,100), bins=150, alpha=0.7, label='DEV', histtype='stepfilled')
plt.hist(r_comp, color='m', range=(0,100), bins=150, alpha=0.7, label= 'COMP', histtype='stepfilled')
plt.hist(r2_rchi2_psf, color='r', range=(0,100), bins=150, alpha=0.9, label='Unmasked CRs', histtype='stepfilled')
plt.title('RChi2 Values (r-filter)')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()

# z_psf histogram
plt.hist(z_psf, color='k', range=(0,100), bins=150, label='PSF', histtype='stepfilled')
plt.hist(z_exp, color='c', range=(0,100), bins=150, alpha=0.7, label='EXP', histtype='stepfilled')
plt.hist(z_dev, color='y', range=(0,100), bins=150, alpha=0.7, label='DEV', histtype='stepfilled')
plt.hist(z_comp, color='m', range=(0,100), bins=150, alpha=0.7, label= 'COMP', histtype='stepfilled')
plt.hist(z2_rchi2_psf, color='r', range=(0,100), bins=200, alpha=0.9, label='Unmasked CRs', histtype='stepfilled')
plt.title('RChi2 Values (z-filter)')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()

i = 0
while i < len(g2_rchi2_psf) and i < len(r2_rchi2_psf) and i < len(z2_rchi2_psf):
    a = g2_rchi2_psf[i]
    b = r2_rchi2_psf[i]
    c = z2_rchi2_psf[i]
    print(a, b, c)
    i += 1
