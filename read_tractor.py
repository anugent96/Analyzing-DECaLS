import sys
image = str(sys.argv[1])
from astropy.io import fits
tractor = fits.open(image)

cols = tractor[1].columns # use to see names of columns
tbl = tractor[1].data # data stored within the table
RA = tbl.field('ra') # RA values
DEC = tbl.field('dec') # DEC values
rE = tbl.field('shapeExp_r') # half-light radius, exponential model
rD = tbl.field('shapeDev_r') # half-light radius, devaucouleurs model
fD = tbl.field('decam_flux') # model flux, decam
eE1 = tbl.field('shapeExp_e1') # ellipticity component 1, exponential model
eD1 = tbl.field('shapeDev_e1') # ellipticity component 1, devaucouleurs model
eE2 = tbl.field('shapeExp_e2') # ellipticity component 2, exponential model
eD2 = tbl.field('shapeDev_e2') # ellipticity component 2, devaucouleurs model
type1 = tbl.field('type') # type: EXP, DEV, PSF, SIMP
Mask = tbl.field('decam_anymask') #check this- any mask (for errors in telescope)
PSF = tbl.field('decam_psfsize') # PSF (Point Spread Function) Size
d_nobs = tbl.field('decam_nobs') # Decam number of exposures



# Necessary imports and definitions
import matplotlib.pyplot as plt
import math
import itertools 
def magnitude(f): # nanomaggies to magnitudes
   if f <= 0:
      m = 35
   else:
      m = 22.5 - 2.5*math.log10(f)
   return m

aMask = [sum(x) for x in Mask]

# Number of observations in z-filter
z_nobs = d_nobs[:,4]
print(min(z_nobs))
print(max(z_nobs))


# RA histogram
plt.hist(RA, histtype='stepfilled', color='c')
plt.title("DECaL RA Histogram")
plt.xlabel("RA Value")
plt.show()

# DEC histogram
plt.hist(DEC, histtype='stepfilled', color='c')
plt.title("DECaL DEC Histogram")
plt.xlabel("DEC Value")
plt.show()

# Radius, Exponential
rE_1 = [word for (word, mask1, mask2) in zip(rE, type1, aMask) if mask1 == 'EXP' and  mask2 == 0]
plt.hist(rE_1, range=(0, 10), histtype='stepfilled', color='c')
plt.title("DECaL Half-Light Radius (Exponential) Histogram")
plt.xlabel("Radius")
plt.yscale('log')
plt.show()

# Radius, de  Vaucouleurs
rD_1 = [word for (word, mask1, mask2) in zip(rD, type1, aMask) if mask1 == 'DEV' and  mask2 == 0]
plt.hist(rD_1, range=(0, 10), histtype='stepfilled', color='c')
plt.title("DECaL Half-Light Radius (de Vaucouleurs) Histogram")
plt.xlabel("Radius")
plt.yscale('log')
plt.show()

# Ellipticity, Exponential
eE1_1 = [word for (word, mask1, mask2) in zip(eE1, type1, aMask) if mask1 == 'EXP' and  mask2 == 0]
plt.hist(eE1_1, bins=60,  histtype='stepfilled', color='c')
plt.title("DECaL Ellipticity Component 1 (Exponential) Histogram")
plt.xlabel("Ellipticity")
plt.show()

eE2_1 = [word for (word, mask1, mask2) in zip(eE2, type1, aMask) if mask1 == 'EXP' and  mask2 == 0]
plt.hist(eE2_1, bins=60, histtype='stepfilled', color='c')
plt.title("DECaL Ellipticity Component 2 (Exponential) Histogram")
plt.xlabel("Ellipticity")
plt.show()
         
# Ellipticity, de Vaucouleurs
eD1_1 = [word for (word, mask1, mask2) in zip(eD1, type1, aMask) if mask1 == 'DEV' and  mask2 == 0]
plt.hist(eD1_1, bins=60,  histtype='stepfilled', color='c')
plt.title("DECaL Ellipticity Component 1 (de Vaucouleurs) Histogram")
plt.xlabel("Ellipticity")
plt.show()

eD2_1 = [word for (word, mask1, mask2) in zip(eD2, type1, aMask) if mask1 == 'DEV' and  mask2 == 0]
plt.hist(eD2_1, bins=60, histtype='stepfilled', color='c')
plt.title("DECaL Ellipticity Component 2 (de Vaucouleurs) Histogram")
plt.xlabel("Ellipticity")
plt.show()

# g-MAG Hist
g = fD[:,1]
gMAG = [magnitude(f) for f in g]
gMAG_1 = [x for (x, mask) in zip(gMAG, aMask) if mask == 0]
plt.hist(gMAG_1, bins=60, histtype='stepfilled', color='c')
plt.title("DeCAL Magnitude Histogram (g-filter)")
plt.xlabel("Magnitude") 
plt.show()

# r-Mag Hist         
r = fD[:,2]
rMAG = [magnitude(f) for f in r]
rMAG_1 = [x for (x, mask) in zip(rMAG, aMask) if mask == 0]
plt.hist(rMAG_1, bins=60, histtype='stepfilled', color='c')
plt.title("DeCAL Magnitude Histogram (r-filter)")
plt.xlabel("Magnitude") 
plt.show()

# z-Mag Hist        
z = fD[:,4]
zMAG = [magnitude(f) for f in z]
zMAG_1 = [x for (x, mask) in zip(zMAG, aMask) if mask == 0]
plt.hist(zMAG_1, bins=60, histtype='stepfilled', color='c')
plt.title("DeCAL Magnitude Histogram (z-filter)")
plt.xlabel("Magnitude") 
plt.show()


