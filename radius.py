import sys
image1 = str(sys.argv[1]) # DECaLS tractor.fits file
image2 = str(sys.argv[2])
image3 = str(sys.argv[3])
image4 = str(sys.argv[4])
image5 = str(sys.argv[5])
image6 = str(sys.argv[6])

from astropy.io import fits
# image 1
tractor1 = fits.open(image1)
tbl1 = tractor1[1].data # data stored within the table
ra1 = tbl1.field('ra') # RA values
dec1 = tbl1.field('dec') # DEC values
rE1 = tbl1.field('shapeExp_r') # half-light radius, exponential model
rD1 = tbl1.field('shapeDev_r') # half-light radius, devaucouleurs model
fD1 = tbl1.field('decam_flux') # model flux, decam
type0 = tbl1.field('type') # type: EXP, DEV, PSF, SIMP
Mask1 = tbl1.field('decam_anymask') #check this- any mask (for errors in telescope)


# image 2
tractor2 = fits.open(image2)
tbl2 = tractor2[1].data # data stored within the table
ra2 = tbl2.field('ra') # RA values
dec2 = tbl2.field('dec') # DEC values
rE2 = tbl2.field('shapeExp_r') # half-light radius, exponential model
rD2 = tbl2.field('shapeDev_r') # half-light radius, devaucouleurs model
fD2 = tbl2.field('decam_flux') # model flux, decam
type2 = tbl2.field('type') # type: EXP, DEV, PSF, SIMP
Mask2 = tbl2.field('decam_anymask') #check this- any mask (for errors in telescope)

# image 3
tractor3 = fits.open(image3)
tbl3 = tractor3[1].data # data stored within the table
ra3 = tbl3.field('ra') # RA values
dec3 = tbl3.field('dec') # DEC values
rE3 = tbl3.field('shapeExp_r') # half-light radius, exponential model
rD3 = tbl3.field('shapeDev_r') # half-light radius, devaucouleurs model
fD3 = tbl3.field('decam_flux') # model flux, decam
type3 = tbl3.field('type') # type: EXP, DEV, PSF, SIMP
Mask3 = tbl3.field('decam_anymask') #check this- any mask (for errors in telescope)

# image 4
tractor4 = fits.open(image4)
tbl4 = tractor4[1].data # data stored within the table
ra4 = tbl4.field('ra') # RA values
dec4 = tbl4.field('dec') # DEC values
rE4 = tbl4.field('shapeExp_r') # half-light radius, exponential model
rD4 = tbl4.field('shapeDev_r') # half-light radius, devaucouleurs model
fD4 = tbl4.field('decam_flux') # model flux, decam
type4 = tbl4.field('type') # type: EXP, DEV, PSF, SIMP
Mask4 = tbl4.field('decam_anymask') #check this- any mask (for errors in telescope)


# image 5
tractor5 = fits.open(image5)
tbl5 = tractor5[1].data # data stored within the table
ra5 = tbl5.field('ra') # RA values
dec5 = tbl5.field('dec') # DEC values
rE5 = tbl5.field('shapeExp_r') # half-light radius, exponential model
rD5 = tbl5.field('shapeDev_r') # half-light radius, devaucouleurs model
fD5 = tbl5.field('decam_flux') # model flux, decam
type5 = tbl5.field('type') # type: EXP, DEV, PSF, SIMP
Mask5 = tbl5.field('decam_anymask') #check this- any mask (for errors in telescope)

# image 6
tractor6 = fits.open(image6)
tbl6 = tractor6[1].data # data stored within the table
ra6 = tbl6.field('ra') # RA values
dec6 = tbl6.field('dec') # DEC values
rE6 = tbl6.field('shapeExp_r') # half-light radius, exponential model
rD6 = tbl6.field('shapeDev_r') # half-light radius, devaucouleurs model
fD6 = tbl6.field('decam_flux') # model flux, decam
type6 = tbl6.field('type') # type: EXP, DEV, PSF, SIMP
Mask6 = tbl6.field('decam_anymask') #check this- any mask (for errors in telescope)


import math
import numpy as np
from numpy import mean
def magnitude(f): # nanomaggies to magnitudes
   if f <= 0:
      m = 35
   else:
      m = 22.5 - 2.5*math.log10(f)
   return m

# Extending lists so that all values are contained in 1 array
import numpy as np
ra = np.append(ra1, ra2)
ra = np.append(ra, ra3)
ra = np.append(ra, ra4)
ra = np.append(ra, ra5)
ra = np.append(ra, ra6)

dec = np.append(dec1, dec2)
dec = np.append(dec, dec3)
dec = np.append(dec, dec4)
dec = np.append(dec, dec5)
dec = np.append(dec, dec6)

rE = np.append(rE1, rE2)
rE = np.append(rE, rE3)
rE = np.append(rE, rE4)
rE = np.append(rE, rE5)
rE = np.append(rE, rE6)

rD = np.append(rD1, rD2)
rD = np.append(rD, rD3)
rD = np.append(rD, rD4)
rD = np.append(rD, rD5)
rD = np.append(rD, rD6)

type1 = np.append(type0, type2)
type1 = np.append(type1, type3)
type1 = np.append(type1, type4)
type1 = np.append(type1, type5)
type1 = np.append(type1, type6)


aMask1 = np.array([sum(x) for x in Mask1]) # sum of each filter mask
aMask2 = np.array([sum(x) for x in Mask2])
aMask3 = np.array([sum(x) for x in Mask3])
aMask4 = np.array([sum(x) for x in Mask4])
aMask5 = np.array([sum(x) for x in Mask5])
aMask6 = np.array([sum(x) for x in Mask6])

aMask = np.append(aMask1, aMask2)
aMask = np.append(aMask, aMask3)
aMask = np.append(aMask, aMask4)
aMask = np.append(aMask, aMask5)
aMask = np.append(aMask, aMask6)

# g, r, and z filters in each of the files
g1 = np.array(fD1[:,1])
r1 = np.array(fD1[:,2])
z1 = np.array(fD1[:,4])

g2 = np.array(fD2[:,1])
r2 = np.array(fD2[:,2])
z2 = np.array(fD2[:,4])

g3 = np.array(fD3[:,1])
r3 = np.array(fD3[:,2])
z3 = np.array(fD3[:,4])

g4 = np.array(fD4[:,1])
r4 = np.array(fD4[:,2])
z4 = np.array(fD4[:,4])

g5 = np.array(fD5[:,1])
r5 = np.array(fD5[:,2])
z5 = np.array(fD5[:,4])

g6 = np.array(fD6[:,1])
r6 = np.array(fD6[:,2])
z6 = np.array(fD6[:,4])

# Extending g-filter list
g = np.append(g1, g2)
g = np.append(g, g3)
g = np.append(g, g4)
g = np.append(g, g5)
g = np.append(g, g6)

# Extending r-filter list
r = np.append(r1, r2)
r = np.append(r, r3)
r = np.append(r, r4)
r = np.append(r, r5)
r = np.append(r, r6)

# Extending z-filter list
z = np.append(z1, z2)
z = np.append(z, z3)
z = np.append(z, z4)
z = np.append(z, z5)
z = np.append(z, z6)

"""
Manipulation of the new lists to get out any maksing errors and convert fluxes to magnitudes.
"""

# Radius (EXP, DEV)
rE_1 = [word for (word, mask1, mask2) in zip(rE, type1, aMask) if (mask1 == 'EXP' or mask1 == 'COMP')and  mask2 == 0]
rD_1 = [word for (word, mask1, mask2) in zip(rD, type1, aMask) if (mask1 == 'DEV' or mask1 == 'COMP') and  mask2 == 0]

# Magnitudes
gMAG = [magnitude(x) for (x, y) in zip(g, aMask) if y == 0]
rMAG = [magnitude(x) for (x, y) in zip(r, aMask) if y == 0]
zMAG = [magnitude(x) for (x, y) in zip(z, aMask) if y == 0]

# Magnitudes at for large radii (>30) (EXP/DEV)
gMagE = [x for (x, y) in zip(gMAG, rE_1) if y > 50]
rMagE = [x for (x, y) in zip(rMAG, rE_1) if y > 50]
zMagE = [x for (x, y) in zip(zMAG, rE_1) if y > 50]

gMagD = [x for (x, y) in zip(gMAG, rD_1) if y > 50]
rMagD = [x for (x, y) in zip(rMAG, rD_1) if y > 50]
zMagD = [x for (x, y) in zip(zMAG, rD_1) if y > 50]

gMagD1 = [x for (x, y) in zip(gMAG, type1) if (y == 'DEV' or y == 'COMP')]
rMagD1 = [x for (x, y) in zip(rMAG, type1) if (y == 'DEV' or y == 'COMP')]
zMagD1 = [x for (x, y) in zip(zMAG, type1) if (y == 'DEV' or y == 'COMP')]

typeD = [x for (x, y) in zip(type1, aMask) if (x == 'DEV' or x == 'COMP') and y == 0]

raC = [x for (x, y) in zip(ra, type1) if y == 'COMP']
print(len(raC))

# Histograms of magnitude at large radii (per type)

import matplotlib.pyplot as plt

plt.hist(gMagE, bins=60, color='k', label='g-filter')
plt.hist(rMagE, bins=60, color='m', alpha=0.7, label='r-filter')
plt.hist(zMagE, bins=60, color='c', alpha=0.4, label='z-filter')
plt.xlabel('Magnitudes')
plt.legend()
plt.title('Magnitude Histogram for Large EXP Radii (r > 50)')
plt.show()

plt.hist(gMagD, bins=60, color='k', label='g-filter')
plt.hist(rMagD, bins=60, color='m', alpha=0.7, label='r-filter')
plt.hist(zMagD, bins=60, color='c', alpha=0.4, label='z-filter')
plt.xlabel('Magnitudes')
plt.legend()
plt.title('Magnitude Histogram for Large DEV Radii (r > 50)')
plt.show()

# RA and DEC printed values for extreme radii
ra_d = [word for (word, mask1, mask2) in zip(ra, type1, aMask) if (mask1 == 'DEV' or mask1 == 'COMP') and  mask2 == 0]
dec_d = [word for (word, mask1, mask2) in zip(dec, type1, aMask) if (mask1 == 'DEV' or mask1 == 'COMP') and  mask2 == 0]
print('Extreme DEV Radii (r > 50): (ra, dec, rDEV)')
i = 0
while i < len(rD_1) and i < len(ra_d) and i < len(dec_d):
    a = rD_1[i]
    b = ra_d[i]
    c = dec_d[i]
    d = gMagD1[i]
    e = rMagD1[i]
    f = zMagD1[i]
    g = typeD[i]
    if a > 50:
        print(b, c, a, d, e, f, g)
    i += 1
