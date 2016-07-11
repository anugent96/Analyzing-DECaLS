import sys
image = str(sys.argv[1])
from astropy.io import fits
tractor = fits.open(image)
tbl = tractor[1].data


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
from numpy import mean
from scipy.stats import moment
def magnitude(f): # nanomaggies to magnitudes
   if f <= 0:
      m = 35
   else:
      m = 22.5 - 2.5*math.log10(f)
   return m

aMask = [sum(x) for x in Mask]

znobs = d_nobs[:,4]

g = fD[:,1]
r = fD[:,2]
z = fD[:,4]

gMag = [magnitude(f) for f in g]
rMag = [magnitude(f) for f in r]
zMag = [magnitude(f) for f in z]

gMag1 = [x for (x,y,w) in zip(gMag, gMag, aMask) if y != 35 and w == 0]
rMag1 = [x for (x,y,w) in zip(rMag, rMag, aMask) if y != 35 and w == 0]
zMag1 = [x for (x,y,w) in zip(zMag, zMag, aMask) if y != 35 and w == 0]

unmasked_obj = [x for (x,y) in zip(z, aMask) if y == 0]
num_unmasked_obj = len(unmasked_obj)

eExp1= [x for (x,y,z) in zip(eE1, aMask, type1) if y == 0 and (z == 'COMP' or z == 'EXP')]
eExp2 = [x for (x,y,z) in zip(eE2, aMask, type1) if y == 0 and (z == 'COMP' or z == 'EXP')]

eDev1= [x for (x,y,z) in zip(eD1, aMask, type1) if y == 0 and (z == 'COMP' or z == 'DEV')]
eDev2 = [x for (x,y,z) in zip(eD2, aMask, type1) if y == 0 and (z == 'COMP' or z == 'DEV')]


print('Average magnitude values (g, r, z):')
print(mean(gMag1), mean(rMag1), mean(zMag1))


print('Maximum magnitude values: ie how faint? (g, r, z)')
print(max(gMag1), max(rMag1), max(zMag1))


print('RA range:')
print(min(RA), max(RA))


print('DEC range:')
print(min(DEC), max(DEC))


print('Exposure range:')
print(min(znobs), max(znobs))


print('PSF Size (g, r, z):')
print(mean(PSF[:,1]), mean(PSF[:,2]), mean(PSF[:,4]))


print('Number of unmasked objects:')
print(num_unmasked_obj)


print('Exponential Ellipticity Skew (avg, moment), (1/2):')
print(mean(eExp1), moment(eExp1, moment=3), mean(eExp2), moment(eExp2, moment=3))


print('de Vaucouleurs Ellipticity Skew (avg, moment), (1/2):')
print(mean(eDev1), moment(eDev1, moment=3), mean(eDev2), moment(eDev2, moment=3))


