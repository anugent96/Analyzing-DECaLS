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

g_rchi2 = rchi2[:,1] # g component of rchi2
r_rchi2 = rchi2[:,2] # r component of rchi2
z_rchi2 = rchi2[:,4] # z component of rchi2

g = fD[:,1] # g flux
r = fD[:,2] # r flux
z = fD[:,4] # z flux

znobs = d_nobs[:,4] # number of observations in z

aMask = [sum(x) for x in mask] # sum of all masking values for 1 object

import matplotlib.pyplot as plt
import math
import itertools

# z rchi2 values in each filter for PSF
g_psf = [x for (x,y,w) in zip(g_rchi2, type1, g_rchi2) if y == 'PSF' and not math.isnan(w)] 
r_psf = [x for (x,y,w) in zip(r_rchi2, type1, r_rchi2) if y == 'PSF' and not math.isnan(w)]
z_psf = [x for (x,y,w) in zip(z_rchi2, type1, z_rchi2) if y == 'PSF' and not math.isnan(w)] 

# z rchi2 values in each filter for EXP
g_exp = [x for (x,y,w) in zip(g_rchi2, type1, g_rchi2) if y == 'EXP' and not math.isnan(w)] 
r_exp = [x for (x,y,w) in zip(r_rchi2, type1, r_rchi2) if y == 'EXP' and not math.isnan(w)]
z_exp = [x for (x,y,w) in zip(z_rchi2, type1, z_rchi2) if y == 'EXP' and not math.isnan(w)] 

# z rchi2 values in each filter for DEV
g_dev = [x for (x,y,w) in zip(g_rchi2, type1, g_rchi2) if y == 'DEV' and not math.isnan(w)] 
r_dev = [x for (x,y,w) in zip(r_rchi2, type1, r_rchi2) if y == 'DEV' and not math.isnan(w)]
z_dev = [x for (x,y,w) in zip(z_rchi2, type1, z_rchi2) if y == 'DEV' and not math.isnan(w)] 

# z rchi2 values in each filter for COMP 
g_comp = [x for (x,y,w) in zip(g_rchi2, type1, g_rchi2) if y == 'COMP' and not math.isnan(w)] 
r_comp = [x for (x,y,w) in zip(r_rchi2, type1, r_rchi2) if y == 'COMP' and not math.isnan(w)]
z_comp = [x for (x,y,w) in zip(z_rchi2, type1, z_rchi2) if y == 'COMP' and not math.isnan(w)] 

# g_psf histogram
plt.hist(g_psf, color='k', range=(0,100), bins=100, label='PSF', histtype='stepfilled')
plt.hist(g_exp, color='c', range=(0,100), bins=100, alpha=0.7, label='EXP', histtype='stepfilled')
plt.hist(g_dev, color='y', range=(0,100), bins=100, alpha=0.7, label='DEV', histtype='stepfilled')
plt.hist(g_comp, color='m', range=(0,100), bins=100, alpha=0.7, label= 'COMP', histtype='stepfilled')
plt.title('RChi2 Values (g-filter)')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()

# r_psf histogram
plt.hist(r_psf, color='k', range=(0,100), bins=100, label='PSF', histtype='stepfilled')
plt.hist(r_exp, color='c', range=(0,100), bins=100, alpha=0.7, label='EXP', histtype='stepfilled')
plt.hist(r_dev, color='y', range=(0,100), bins=100, alpha=0.7, label='DEV', histtype='stepfilled')
plt.hist(r_comp, color='m', range=(0,100), bins=100, alpha=0.7, label= 'COMP', histtype='stepfilled')
plt.title('RChi2 Values (r-filter)')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()

# z_psf histogram
plt.hist(z_psf, color='k', range=(0,100), bins=100, label='PSF', histtype='stepfilled')
plt.hist(z_exp, color='c', range=(0,100), bins=100, alpha=0.7, label='EXP', histtype='stepfilled')
plt.hist(z_dev, color='y', range=(0,100), bins=100, alpha=0.7, label='DEV', histtype='stepfilled')
plt.hist(z_comp, color='m', range=(0,100), bins=100, alpha=0.7, label= 'COMP', histtype='stepfilled')
plt.title('RChi2 Values (z-filter)')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()