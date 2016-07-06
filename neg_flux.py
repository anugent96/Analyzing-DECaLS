"""
Finds RA, dec, fluxes (g,r,z,w1), type, masking value, and z-n-obs for when z-flux is positive and all other fluxes are negative. Works
with any 1 fits file.
"""

import sys
image1 = str(sys.argv[1])
from astropy.io import fits
tractor = fits.open(image1)

tbl = tractor[1].data 
RA = tbl.field('ra') # RA values
DEC = tbl.field('dec') # DEC values
fW = tbl.field('wise_flux') # model flux, WISE
fD = tbl.field('decam_flux') # model flux, decam
Mask = tbl.field('decam_anymask') #check this- any mask (for errors in telescope)
d_nobs = tbl.field('decam_nobs') # decam nobs
Type = tbl.field('type') # type of object

# Necessary Imports/ Functions
g = fD[:,1] # g-flux
r = fD[:,2] # r-flux
z = fD[:,4] # z-flux
w1 = fW[:,0] # w1 (WISE) flux
z_nobs = d_nobs[:,4] # number of observations with z-filter
aMask = [sum(x) for x in Mask] # see if any individual filter error is there for a candidate

# To print negative fluxes
RA1 = [x for (x, mask1, mask2) in zip(RA, w1, z) if mask1 < 0 and mask2 > 0] # when z-flux is positive and w1-flux negative
DEC1 = [x for (x, mask1, mask2) in zip(DEC, w1, z) if mask1 < 0 and mask2 > 0]
g1 = [x for (x, mask1, mask2) in zip(g, w1, z) if mask1 < 0 and mask2 > 0]
r1 = [x for (x, mask1, mask2) in zip(r, w1, z) if mask1 < 0 and mask2 > 0]
z1 = [x for (x, mask1, mask2) in zip(z, w1, z) if mask1 < 0 and mask2 > 0]
w1_1 = [x for (x, mask1, mask2) in zip(w1, w1, z) if mask1 < 0 and mask2 > 0]
aMask1 = [x for (x, mask1, mask2) in zip(aMask, w1, z) if mask1 < 0 and mask2 > 0]
type1 = [x for (x, mask1, mask2) in zip (Type, w1, z ) if mask1 < 0 and mask2 > 0]
z_nobs1 = [x for (x, mask1, mask2) in zip (z_nobs, w1, z ) if mask1 < 0 and mask2 > 0]

RA2 = [x for (x, mask1, mask2) in zip(RA1, g1, r1) if mask1 < 0 and mask2 < 0] # when above conditions are met and g and r fluxes are negative
DEC2 = [x for (x, mask1, mask2) in zip(DEC1, g1, r1) if mask1 < 0 and mask2 < 0]
g2 = [x for (x, mask1, mask2) in zip(g1, g1, r1) if mask1 < 0 and mask2 < 0]
r2 = [x for (x, mask1, mask2) in zip(r1, g1, r1) if mask1 < 0 and mask2 < 0]
z2 = [x for (x, mask1, mask2) in zip(z1, g1, r1) if mask1 < 0 and mask2 < 0]
w1_2 = [x for (x, mask1, mask2) in zip(w1_1, g1, r1) if mask1 < 0 and mask2 < 0]
aMask2 = [x for (x, mask1, mask2) in zip(aMask1, g1, r1) if mask1 < 0 and mask2 < 0]
type2 = [x for (x, mask1, mask2) in zip (type1, g1, r1) if mask1 < 0 and mask2 < 0]
z_nobs2 = [x for (x, mask1, mask2) in zip (z_nobs1, g1, r1) if mask1 < 0 and mask2 < 0]

g_neg = [x for x in g if x < 0]
r_neg = [x for x in r if x < 0]
z_neg = [word for word in z if word < 0]
print(len(g_neg), len(r_neg)) # Prints TOTAL number of negative fluxes in g, r, and z filters
print(len(z_neg))

"""
To print all RA, dec, fluxes, masking value, type, and z-n-obs for object that match the conditions (z pos and everything else neg)
"""

print('RA , DEC, g, r, z, w1, mask, type, z-nobs')
i = 0
while i < len(RA2) and i < len(DEC2) and i < len(g2) and i < len(r2) and i < len(z2) and i < len(w1_2) and i < len(aMask2):
    x = RA2[i]
    y = DEC2[i]
    z = g2[i]
    w = r2[i]
    a = z2[i]
    b = w1_2[i]
    c = aMask2[i]
    d = type2[i]
    f = z_nobs2[i]
    print(x, y, z, w, a, b, c, d, f)
    i +=1

