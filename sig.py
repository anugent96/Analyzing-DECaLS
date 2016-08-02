""" 
Takes any 1 fits fle and prints RA, dec, fluxes (g, r, z, w1), masking value, type, number of observations in z, and rchi2 value for when z flux
is significantly positive (>2) and all other fluxes are not significant (<1).
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
d_nobs = tbl.field('decam_nobs') # number of observations taken in each filter
Type = tbl.field('type') # type of object
ivd = tbl.field('decam_flux_ivar') # inverse variance for decam (1/sigma squared)
ivw = tbl.field('wise_flux_ivar')  # inverse variance for wise
rchi2 = tbl.field('decam_rchi2') # decam rchi2 value

# Necessary Imports/ Functions
import numpy as np
g = np.array(fD[:,1]) # g-flux
r = np.array(fD[:,2]) # r-flux
z = np.array(fD[:,4]) # z-flux
w1 = np.array(fW[:,0]) # w1-flux
z_nobs = d_nobs[:,4] # number of observations in z
g_sig2 = np.array(ivd[:,1]) # g inverse variance
r_sig2 = np.array(ivd[:,2]) # r inverse variance
z_sig2 = np.array(ivd[:,4]) # z inverse variance
w1_sig2 = np.array(ivw[:,0]) # w1 inverse variance
aMask = [sum(x) for x in Mask] # see if any individual filter error is there for a candidate
z_rchi2 = rchi2[:,4] # rchi2 value in z filter


"""
Flux times the squareroot of its corresponding inverse variance. Gives a better value for flux that can be used to find fluxes that
are significant and not significant (ie: objects that exist and don't exist)
"""
g_test = g*(g_sig2**0.5)  
r_test = r*(r_sig2**0.5)
z_test = z*(z_sig2**0.5)
w1_test = w1*(w1_sig2**0.5)


# To print negative fluxes
"""
Finds RA, DEC, fluxes, masking values, type, rchi2 for when object is significantly visible in only 1 image.
"""

# Possible CRs in z
RA1 = [x for (x, mask1, mask2) in zip(RA, w1_test, z_test) if mask1 < 1 and mask2 > 2]
DEC1 = [x for (x, mask1, mask2) in zip(DEC, w1_test, z_test) if mask1 < 1 and mask2 > 2]
g1 = [x for (x, mask1, mask2) in zip(g_test, w1_test, z_test) if mask1 < 1 and mask2 > 2]
r1 = [x for (x, mask1, mask2) in zip(r_test, w1_test, z_test) if mask1 < 1 and mask2 > 2]
z1 = [x for (x, mask1, mask2) in zip(z_test, w1_test, z_test) if mask1 < 1 and mask2 > 2]
w1_1 = [x for (x, mask1, mask2) in zip(w1_test, w1_test, z_test) if mask1 < 1 and mask2 > 2]
aMask1 = [x for (x, mask1, mask2) in zip(aMask, w1_test, z_test) if mask1 < 1 and mask2 > 2]
type1 = [x for (x, mask1, mask2) in zip (Type, w1_test, z_test) if mask1 < 1 and mask2 > 2]
z_nobs1 = [x for (x, mask1, mask2) in zip (z_nobs, w1_test, z_test) if mask1 < 1 and mask2 > 2]
z1_rchi2 = [x for (x, mask1, mask2) in zip (z_rchi2, w1_test, z_test) if mask1 < 1 and mask2 > 2]
      
RA2 = [x for (x, mask1, mask2) in zip(RA1, g1, r1) if mask1 < 1 and mask2 < 1]
DEC2 = [x for (x, mask1, mask2) in zip(DEC1, g1, r1) if mask1 < 1 and mask2 < 1]
g2 = [x for (x, mask1, mask2) in zip(g1, g1, r1) if mask1 < 1 and mask2 < 1]
r2 = [x for (x, mask1, mask2) in zip(r1, g1, r1) if mask1 < 1 and mask2 < 1]
z2 = [x for (x, mask1, mask2) in zip(z1, g1, r1) if mask1 < 1 and mask2 < 1]
w1_2 = [x for (x, mask1, mask2) in zip(w1_1, g1, r1) if mask1 < 1 and mask2 < 1]
aMask2 = [x for (x, mask1, mask2) in zip(aMask1, g1, r1) if mask1 < 1 and mask2 < 1]
type2 = [x for (x, mask1, mask2) in zip (type1, g1, r1) if mask1 < 1 and mask2 < 1]
z_nobs2 = [x for (x, mask1, mask2) in zip (z_nobs1, g1, r1) if mask1 < 1 and mask2 < 1]
z2_rchi2 = [x for (x, mask1, mask2) in zip (z1_rchi2, g1, r1) if mask1 < 1 and mask2 < 1]

# Possible CRs in g
RA3 = [x for (x, mask1, mask2) in zip(RA, w1_test, g_test) if mask1 < 1 and mask2 > 2]
aMask3 = [x for (x, mask1, mask2) in zip(aMask, w1_test, g_test) if mask1 < 1 and mask2 > 2]
g3 = [x for (x, mask1, mask2) in zip(g_test, w1_test, g_test) if mask1 < 1 and mask2 > 2]
r3 = [x for (x, mask1, mask2) in zip(r_test, w1_test, g_test) if mask1 < 1 and mask2 > 2]
z3 = [x for (x, mask1, mask2) in zip(z_test, w1_test, g_test) if mask1 < 1 and mask2 > 2]
w1_3 = [x for (x, mask1, mask2) in zip(w1_test, w1_test, g_test) if mask1 < 1 and mask2 > 2]

      
RA4 = [x for (x, mask1, mask2) in zip(RA3, z3, r3) if mask1 < 1 and mask2 < 1]
aMask4 = [x for (x, mask1, mask2) in zip(aMask3, z3, r3) if mask1 < 1 and mask2 < 1]
RA5 = [x for (x, y) in zip(RA4, aMask4) if y == 0]

# Possible CRs in r
RA6 = [x for (x, mask1, mask2) in zip(RA, w1_test, r_test) if mask1 < 1 and mask2 > 2]
aMask6 = [x for (x, mask1, mask2) in zip(aMask, w1_test, r_test) if mask1 < 1 and mask2 > 2]
g6 = [x for (x, mask1, mask2) in zip(g_test, w1_test, r_test) if mask1 < 1 and mask2 > 2]
r6 = [x for (x, mask1, mask2) in zip(r_test, w1_test, r_test) if mask1 < 1 and mask2 > 2]
z6 = [x for (x, mask1, mask2) in zip(z_test, w1_test, r_test) if mask1 < 1 and mask2 > 2]
w1_6 = [x for (x, mask1, mask2) in zip(w1_test, w1_test, r_test) if mask1 < 1 and mask2 > 2]

      
RA7 = [x for (x, mask1, mask2) in zip(RA6, z6, g6) if mask1 < 1 and mask2 < 1]
aMask7 = [x for (x, mask1, mask2) in zip(aMask6, z6, g6) if mask1 < 1 and mask2 < 1]
RA8 = [x for (x, y) in zip(RA7, aMask7) if y == 0]


print('RA , DEC, g, r, z, w1, mask, type, z-nobs, z-rchi2')
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
    g = z2_rchi2[i]
    if c == 0:
        print(x, y, z, w, a, b, c, d, f, g)    
    i +=1

 
print('Num total unmasked objects detected:')
RA_unmask = [x for (x,y) in zip(RA, aMask) if y == 0]
print(len(RA_unmask))

print('Num possible unmasked CRs in z=')
print(len(RA2))

print('Num possible unmasked CRs in g=')
print(len(RA5))

print('Num possible unmasked CRs in r=')
print(len(RA8))


zPsf1 = [x for (x, y, z) in zip(z2_rchi2, aMask, z2_rchi2) if y == 0 and z < 1]
print('Number of possible CRs with z-RChi2 value under 1:')
print(len(zPsf1))

