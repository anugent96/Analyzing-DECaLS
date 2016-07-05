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
d_nobs = tbl.field('decam_nobs')
Type = tbl.field('type') # type of object
ivd = tbl.field('decam_flux_ivar') # inverse variance (1/sigma squared)
ivw = tbl.field('wise_flux_ivar')

# Necessary Imports/ Functions
import numpy as np
g = np.array(fD[:,1])
r = np.array(fD[:,2])
z = np.array(fD[:,4])
w1 = np.array(fW[:,0])
z_nobs = d_nobs[:,4]
g_sig2 = np.array(ivd[:,1])
r_sig2 = np.array(ivd[:,2])
z_sig2 = np.array(ivd[:,4])
w1_sig2 = np.array(ivw[:,0])
aMask = [sum(x) for x in Mask] # see if any individual filter error is there for a candidate

g_test = g*(g_sig2**0.5) 
r_test = r*(r_sig2**0.5)
z_test = z*(z_sig2**0.5)
w1_test = w1*(w1_sig2**0.5)


# To print negative fluxes
RA1 = [x for (x, mask1, mask2) in zip(RA, w1_test, z_test) if mask1 < 1 and mask2 > 2]
DEC1 = [x for (x, mask1, mask2) in zip(DEC, w1_test, z_test) if mask1 < 1 and mask2 > 2]
g1 = [x for (x, mask1, mask2) in zip(g_test, w1_test, z_test) if mask1 < 1 and mask2 > 2]
r1 = [x for (x, mask1, mask2) in zip(r_test, w1_test, z_test) if mask1 < 1 and mask2 > 2]
z1 = [x for (x, mask1, mask2) in zip(z_test, w1_test, z_test) if mask1 < 1 and mask2 > 2]
w1_1 = [x for (x, mask1, mask2) in zip(w1_test, w1_test, z_test) if mask1 < 1 and mask2 > 2]
aMask1 = [x for (x, mask1, mask2) in zip(aMask, w1_test, z_test) if mask1 < 1 and mask2 > 2]
type1 = [x for (x, mask1, mask2) in zip (Type, w1_test, z_test) if mask1 < 1 and mask2 > 2]
z_nobs1 = [x for (x, mask1, mask2) in zip (z_nobs, w1_test, z_test) if mask1 < 1 and mask2 > 2]

      
RA2 = [x for (x, mask1, mask2) in zip(RA1, g1, r1) if mask1 < 1 and mask2 < 1]
DEC2 = [x for (x, mask1, mask2) in zip(DEC1, g1, r1) if mask1 < 1 and mask2 < 1]
g2 = [x for (x, mask1, mask2) in zip(g1, g1, r1) if mask1 < 1 and mask2 < 1]
r2 = [x for (x, mask1, mask2) in zip(r1, g1, r1) if mask1 < 1 and mask2 < 1]
z2 = [x for (x, mask1, mask2) in zip(z1, g1, r1) if mask1 < 1 and mask2 < 1]
w1_2 = [x for (x, mask1, mask2) in zip(w1_1, g1, r1) if mask1 < 1 and mask2 < 1]
aMask2 = [x for (x, mask1, mask2) in zip(aMask1, g1, r1) if mask1 < 1 and mask2 < 1]
type2 = [x for (x, mask1, mask2) in zip (type1, g1, r1) if mask1 < 1 and mask2 < 1]
z_nobs2 = [x for (x, mask1, mask2) in zip (z_nobs1, g1, r1) if mask1 < 1 and mask2 < 1]


# /Users/anuge96/Downloads/ds9.darwinelcapitan.7.4/ds9 -zscale ~/Desktop/Fits\ Images/decals-1791p015-image-z.fits  -region ds9.reg    >> to circle all of these points
# Region file format: DS9 version 4.1
# global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
# fk5

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
    if c == 0:
        print(x, y, z, w, a, b, c, d, f)    
    i +=1

    
    
print('Num possible unmasked CRs=', len(RA2))    


