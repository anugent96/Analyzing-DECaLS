"""
Finds "extreme" (e < -.75 and e > .75) ellipticities and plots histograms for their radii and magnitudes. Used to see if extreme 
ellipticities correspond to a specific type of galaxy or star (eg: small, faint)
"""


import sys
image1 = str(sys.argv[1])
image2 = str(sys.argv[2])
from astropy.io import fits
tractor = fits.open(image1)
hdulist = fits.open(image2)
cols = tractor[1].columns # use to see names of columns
tbl = tractor[1].data # data stored within the table
RA_A = tbl.field('ra') # RA values
DEC_A = tbl.field('dec') # DEC values
rE_A = tbl.field('shapeExp_r') # half-light radius, exponential model
rD_A = tbl.field('shapeDev_r') # half-light radius, devaucouleurs model
fD_A = tbl.field('decam_flux') # model flux, decam
eE1_A = tbl.field('shapeExp_e1') # ellipticity component 1, exponential model
eD1_A = tbl.field('shapeDev_e1') # ellipticity component 1, devaucouleurs model
eE2_A = tbl.field('shapeExp_e2') # ellipticity component 2, exponential model
eD2_A = tbl.field('shapeDev_e2') # ellipticity component 2, devaucouleurs model
type_A = tbl.field('type') # type: EXP, DEV, PSF, SIMP
Mask_A = tbl.field('decam_anymask') #check this- any mask (for errors in telescope)
PSF_A = tbl.field('decam_psfsize') # PSF (Point Spread Function) Size

columns = hdulist[1].columns # use to see names of columns
data = hdulist[1].data # data stored within the table
RA_B = data.field('ra') # RA values
DEC_B = data.field('dec') # DEC values
rE_B = data.field('shapeExp_r') # half-light radius, exponential model
rD_B = data.field('shapeDev_r') # half-light radius, devaucouleurs model
fD_B = data.field('decam_flux') # model flux, decam
eE1_B = data.field('shapeExp_e1') # ellipticity component 1, exponential model
eD1_B = data.field('shapeDev_e1') # ellipticity component 1, devaucouleurs model
eE2_B = data.field('shapeExp_e2') # ellipticity component 2, exponential model
eD2_B = data.field('shapeDev_e2') # ellipticity component 2, devaucouleurs model
type_B = data.field('type') # type: EXP, DEV, PSF, SIMP
Mask_B = data.field('decam_anymask') #check this- any mask (for errors in telescope)
PSF_B = data.field('decam_psfsize') # PSF (Point Spread Function) Size

import matplotlib.pyplot as plt
import math
import itertools 
def magnitude(f): # nanomaggies to magnitudes
   if f <= 0: # if the flux is less then 0, magnitude value will be set to 35 so that it is separated from the rest of the data
      m = 35
   else:
      m = 22.5 - 2.5*math.log10(f)
   return m

aMask_A = [sum(x) for x in Mask_A] # see if any individual filter error is there for a candidate
aMask_B = [sum(x) for x in Mask_B]

"""
DEV and EXP galaxies' ellipticity and radii will be calculated separately as to not mix-up the data. We also only
want to look unbiased data, so all objects with a masking value > 0 will be thrown out
"""

rE_A1 = [word for (word, mask1, mask2) in zip(rE_A, type_A, aMask_A) if mask1 == 'EXP' and  mask2 == 0]
rE_B1 = [word for (word, mask1, mask2) in zip(rE_B, type_B, aMask_B) if mask1 == 'EXP' and  mask2 == 0]
rD_A1 = [word for (word, mask1, mask2) in zip(rD_A, type_A, aMask_A) if mask1 == 'DEV' and  mask2 == 0]
rD_B1 = [word for (word, mask1, mask2) in zip(rD_B, type_B, aMask_B) if mask1 == 'DEV' and  mask2 == 0]
eE1_A1 = [word for (word, mask1, mask2) in zip(eE1_A, type_A, aMask_A) if mask1 == 'EXP' and  mask2 == 0]
eE1_B1 = [word for (word, mask1, mask2) in zip(eE1_B, type_B, aMask_B) if mask1 == 'EXP' and  mask2 == 0]
eE2_A2 = [word for (word, mask1, mask2) in zip(eE2_A, type_A, aMask_A) if mask1 == 'EXP' and  mask2 == 0]
eE2_B2 = [word for (word, mask1, mask2) in zip(eE2_B, type_B, aMask_B) if mask1 == 'EXP' and  mask2 == 0]
eD1_A1 = [word for (word, mask1, mask2) in zip(eD1_A, type_A, aMask_A) if mask1 == 'DEV' and  mask2 == 0]
eD1_B1 = [word for (word, mask1, mask2) in zip(eD1_B, type_B, aMask_B) if mask1 == 'DEV' and  mask2 == 0]
eD2_A2 = [word for (word, mask1, mask2) in zip(eD2_A, type_A, aMask_A) if mask1 == 'DEV' and  mask2 == 0]
eD2_B2 = [word for (word, mask1, mask2) in zip(eD2_B, type_B, aMask_B) if mask1 == 'DEV' and  mask2 == 0]

"""
Finding g, r, and z fluxes and their respective magnitudes.
"""

g_A = fD_A[:,1]
g_B = fD_B[:,1]
gMAG_A = [magnitude(f) for f in g_A]
gMAG_B = [magnitude(f) for f in g_B]
gMAG_A1 = [x for (x, mask) in zip(gMAG_A, aMask_A) if mask == 0]
gMAG_B1 = [x for (x, mask) in zip(gMAG_B, aMask_B) if mask == 0]
r_A = fD_A[:,2]
r_B = fD_B[:,2]
rMAG_A = [magnitude(f) for f in r_A]
rMAG_B = [magnitude(f) for f in r_B]
rMAG_A1 = [x for (x, mask) in zip(rMAG_A, aMask_A) if mask == 0]
rMAG_B1 = [x for (x, mask) in zip(rMAG_B, aMask_B) if mask == 0]
z_A = fD_A[:,4]
z_B = fD_B[:,4]
zMAG_A = [magnitude(f) for f in z_A]
zMAG_B = [magnitude(f) for f in z_B]
zMAG_A1 = [x for (x, mask) in zip(zMAG_A, aMask_A) if mask == 0]
zMAG_B1 = [x for (x, mask) in zip(zMAG_B, aMask_B) if mask == 0]

#Exponential Graphs
rE1_A_ellipse = [x for (x, mask1) in zip(rE_A1, eE1_A1) if mask1 >= 0.75 or mask1<= -0.75]  # conditions for extreme ellipticities
rE1_B_ellipse = [x for (x, mask1) in zip(rE_B1, eE1_B1) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(rE1_A_ellipse, histtype='stepfilled', color='r')
plt.hist(rE1_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity Radius Histogram (Exponential 1)")
plt.xlabel("Radius (e>0.75, e<-0.75)")
plt.show()

rE2_A_ellipse = [x for (x, mask1) in zip(rE_A1, eE2_A2) if mask1 >= 0.75 or mask1<= -0.75] 
rE2_B_ellipse = [x for (x, mask1) in zip(rE_B1, eE2_B2) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(rE2_A_ellipse, histtype='stepfilled', color='r')
plt.hist(rE2_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity Radius Histogram (Exponential 2)")
plt.xlabel("Radius (e>0.75, e<-0.75)")
plt.show()

g1_A_ellipse = [x for (x, mask1) in zip(gMAG_A1, eE1_A1) if mask1 >= 0.75 or mask1<= -0.75] 
g1_B_ellipse = [x for (x, mask1) in zip(gMAG_B1, eE1_B1) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(g1_A_ellipse, histtype='stepfilled', color='r')
plt.hist(g1_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity g-Magnitude Histogram (Exponential 1)")
plt.xlabel("Magnitude (e>0.75, e<-0.75)")
plt.show()

g2_A_ellipse = [x for (x, mask1) in zip(gMAG_A1, eE2_A2) if mask1 >= 0.75 or mask1<= -0.75] 
g2_B_ellipse = [x for (x, mask1) in zip(gMAG_B1, eE2_B2) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(g2_A_ellipse, histtype='stepfilled', color='r')
plt.hist(g2_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity g-Magnitude Histogram (Exponential 2)")
plt.xlabel("Magnitude (e>0.75, e<-0.75)")
plt.show()

r1_A_ellipse = [x for (x, mask1) in zip(rMAG_A1, eE1_A1) if mask1 >= 0.75 or mask1<= -0.75] 
r1_B_ellipse = [x for (x, mask1) in zip(rMAG_B1, eE1_B1) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(r1_A_ellipse, histtype='stepfilled', color='r')
plt.hist(r1_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity r-Magnitude Histogram (Exponential 1)")
plt.xlabel("Magnitude (e>0.75, e<-0.75)")
plt.show()

r2_A_ellipse = [x for (x, mask1) in zip(rMAG_A1, eE2_A2) if mask1 >= 0.75 or mask1<= -0.75] 
r2_B_ellipse = [x for (x, mask1) in zip(rMAG_B1, eE2_B2) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(r2_A_ellipse, histtype='stepfilled', color='r')
plt.hist(r2_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity r-Magnitude Histogram (Exponential 2)")
plt.xlabel("Magnitude (e>0.75, e<-0.75)")
plt.show()

z1_A_ellipse = [x for (x, mask1) in zip(zMAG_A1, eE1_A1) if mask1 >= 0.75 or mask1<= -0.75] 
z1_B_ellipse = [x for (x, mask1) in zip(zMAG_B1, eE1_B1) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(g1_A_ellipse, histtype='stepfilled', color='r')
plt.hist(g1_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity z-Magnitude Histogram (Exponential 1)")
plt.xlabel("Magnitude (e>0.75, e<-0.75)")
plt.show()

z2_A_ellipse = [x for (x, mask1) in zip(zMAG_A1, eE2_A2) if mask1 >= 0.75 or mask1<= -0.75] 
z2_B_ellipse = [x for (x, mask1) in zip(zMAG_B1, eE2_B2) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(z2_A_ellipse, histtype='stepfilled', color='r')
plt.hist(z2_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity z-Magnitude Histogram (Exponential 2)")
plt.xlabel("Magnitude (e>0.75, e<-0.75)")
plt.show()

# de Vaucouleurs
rD1_A_ellipse = [x for (x, mask1) in zip(rD_A1, eD1_A1) if mask1 >= 0.75 or mask1<= -0.75] 
rD1_B_ellipse = [x for (x, mask1) in zip(rD_B1, eD1_B1) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(rD1_A_ellipse, histtype='stepfilled', color='r')
plt.hist(rD1_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity Radius Histogram (de Vaucouleurs 1)")
plt.xlabel("Radius (e>0.75, e<-0.75)")
plt.show()

rD2_A_ellipse = [x for (x, mask1) in zip(rD_A1, eD2_A2) if mask1 >= 0.75 or mask1<= -0.75] 
rD2_B_ellipse = [x for (x, mask1) in zip(rD_B1, eD2_B2) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(rD2_A_ellipse, histtype='stepfilled', color='r')
plt.hist(rD2_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity Radius Histogram (de Vaucouleurs 2)")
plt.xlabel("Radius (e>0.75, e<-0.75)")
plt.show()

g1_A_ellipse = [x for (x, mask1) in zip(gMAG_A1, eD1_A1) if mask1 >= 0.75 or mask1<= -0.75] 
g1_B_ellipse = [x for (x, mask1) in zip(gMAG_B1, eD1_B1) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(g1_A_ellipse, histtype='stepfilled', color='r')
plt.hist(g1_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity g-Magnitude Histogram (de Vaucouleurs 1)")
plt.xlabel("Magnitude (e>0.75, e<-0.75)")
plt.show()

g2_A_ellipse = [x for (x, mask1) in zip(gMAG_A1, eD2_A2) if mask1 >= 0.75 or mask1<= -0.75] 
g2_B_ellipse = [x for (x, mask1) in zip(gMAG_B1, eD2_B2) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(g2_A_ellipse, histtype='stepfilled', color='r')
plt.hist(g2_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity g-Magnitude Histogram (de Vaucouleurs 2)")
plt.xlabel("Magnitude (e>0.75, e<-0.75)")
plt.show()

r1_A_ellipse = [x for (x, mask1) in zip(rMAG_A1, eD1_A1) if mask1 >= 0.75 or mask1<= -0.75] 
r1_B_ellipse = [x for (x, mask1) in zip(rMAG_B1, eD1_B1) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(r1_A_ellipse, histtype='stepfilled', color='r')
plt.hist(r1_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity r-Magnitude Histogram (de Vaucouleurs 1)")
plt.xlabel("Magnitude (e>0.75, e<-0.75)")
plt.show()

r2_A_ellipse = [x for (x, mask1) in zip(rMAG_A1, eD2_A2) if mask1 >= 0.75 or mask1<= -0.75] 
r2_B_ellipse = [x for (x, mask1) in zip(rMAG_B1, eD2_B2) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(r2_A_ellipse, histtype='stepfilled', color='r')
plt.hist(r2_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity r-Magnitude Histogram (de Vaucouleurs 2)")
plt.xlabel("Magnitude (e>0.75, e<-0.75)")
plt.show()

z1_A_ellipse = [x for (x, mask1) in zip(zMAG_A1, eD1_A1) if mask1 >= 0.75 or mask1<= -0.75] 
z1_B_ellipse = [x for (x, mask1) in zip(zMAG_B1, eD1_B1) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(g1_A_ellipse, histtype='stepfilled', color='r')
plt.hist(g1_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity z-Magnitude Histogram (de Vaucouleurs 1)")
plt.xlabel("Magnitude (e>0.75, e<-0.75)")
plt.show()

z2_A_ellipse = [x for (x, mask1) in zip(zMAG_A1, eD2_A2) if mask1 >= 0.75 or mask1<= -0.75] 
z2_B_ellipse = [x for (x, mask1) in zip(zMAG_B1, eD2_B2) if mask1 >= 0.75 or mask1<= -0.75]
plt.hist(z2_A_ellipse, histtype='stepfilled', color='r')
plt.hist(z2_B_ellipse, histtype='stepfilled', color='b', alpha=0.7)
plt.title("Extreme Ellipticity z-Magnitude Histogram (de Vaucouleurs 2)")
plt.xlabel("Magnitude (e>0.75, e<-0.75)")
plt.show()
