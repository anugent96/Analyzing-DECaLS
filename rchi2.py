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

g_psf = [x for (x,y,w) in zip(g_rchi2, type1, g_rchi2) if y == 'PSF' and not math.isnan(w)] # g rchi2 values for PSF types
r_psf = [x for (x,y,w) in zip(r_rchi2, type1, r_rchi2) if y == 'PSF' and not math.isnan(w)]# r rchi2 values for PSF types
z_psf = [x for (x,y,w) in zip(z_rchi2, type1, z_rchi2) if y == 'PSF' and not math.isnan(w)] # z rchi2 values for PSF types

RA_psf = [x for (x,y,w) in zip(RA, type1, z_rchi2) if y == 'PSF' and not math.isnan(w)]
dec_psf = [x for (x,y,w) in zip(DEC, type1, z_rchi2) if y == 'PSF' and not math.isnan(w)]
aMask_psf = [x for (x,y,w) in zip(aMask, type1, z_rchi2) if y == 'PSF' and not math.isnan(w)]

"""
# g_psf histogram
plt.hist(g_psf, color='c', range=(0,100), bins=100)
plt.title('RChi2 Values for PSF Objects (g-filter)')
plt.xscale('log')
plt.yscale('log')
plt.show()


# r_psf histogram
plt.hist(r_psf, color='c', bins=100, range=(0,100))
plt.title('RChi2 Values for PSF Objects (r-filter)')
plt.xscale('log')
plt.yscale('log')
plt.show()

# z_psf histogram
plt.hist(z_psf, color='c', bins=100, range=(0,100))
plt.title('RChi2 Values for PSF Objects (z-filter)')
plt.xscale('log')
plt.yscale('log')
plt.show()
"""
print('length of list:', len(RA))
print('number of PSFs:', len(z_psf))

print('z-rChi2, RA, DEC, masking value:')

i = 0
while i < len(z_psf) and i < len(RA_psf) and len(dec_psf) and len(aMask_psf):
    a = z_psf[i]
    b = RA_psf[i]
    c = dec_psf[i]
    d = aMask_psf[i]
    if a < 1 and d == 0:
        print(a, b, c, d)
    i += 1   

zPsf1 = [x for (x, y, z) in zip(z_psf, aMask, z_psf) if y == 0 and z < 1]
print('Number of objects with RChi2 value under 1:')
print(len(zPsf1))
    