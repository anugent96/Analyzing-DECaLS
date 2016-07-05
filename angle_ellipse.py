import sys
image1 = str(sys.argv[1])
image2 = str(sys.argv[2])
from astropy.io import fits
tractor1 = fits.open(image1)
tractor2 = fits.open(image2)

tbl = tractor1[1].data 
data = tractor2[1].data

eE1_A = tbl.field('shapeExp_e1') # ellipticity component 1, exponential model
eD1_A = tbl.field('shapeDev_e1') # ellipticity component 1, devaucouleurs model
eE2_A = tbl.field('shapeExp_e2') # ellipticity component 2, exponential model
eD2_A = tbl.field('shapeDev_e2') # ellipticity component 2, devaucouleurs model
Mask_A = tbl.field('decam_anymask') #check this- any mask (for errors in telescope)

eE1_B = data.field('shapeExp_e1') # ellipticity component 1, exponential model
eD1_B = data.field('shapeDev_e1') # ellipticity component 1, devaucouleurs model
eE2_B = data.field('shapeExp_e2') # ellipticity component 2, exponential model
eD2_B = data.field('shapeDev_e2') # ellipticity component 2, devaucouleurs model
Mask_B = data.field('decam_anymask') #check this- any mask (for errors in telescope)


import matplotlib.pyplot as plt
from math import atan
from math import degrees
import itertools
import numpy
from scipy.stats import moment
from operator import truediv
aMask_A = [sum(x) for x in Mask_A] # see if any individual filter error is there for a candidate
aMask_B = [sum(x) for x in Mask_B]

def orientation(x):
    angle = 0.5*atan(x)
    return degrees(angle)

eE1_A1 = [x for (x, mask1, mask2) in zip(eE1_A, eE1_A, aMask_A) if mask1 != 0 and mask2 == 0]
eE2_A1 = [x for (x, mask1, mask2) in zip(eE2_A, eE1_A, aMask_A) if mask1 != 0 and mask2 == 0]
eD1_A1 = [x for (x, mask1, mask2) in zip(eD1_A, eD1_A, aMask_A) if mask1 != 0 and mask2 == 0]
eD2_A1 = [x for (x, mask1, mask2) in zip(eD2_A, eD1_A, aMask_A) if mask1 != 0 and mask2 == 0]

eE1_B1 = [x for (x, mask1, mask2) in zip(eE1_B, eE1_B, aMask_B) if mask1 != 0 and mask2 == 0]
eE2_B1 = [x for (x, mask1, mask2) in zip(eE2_B, eE1_B, aMask_B) if mask1 != 0 and mask2 == 0]
eD1_B1 = [x for (x, mask1, mask2) in zip(eD1_B, eD1_B, aMask_B) if mask1 != 0 and mask2 == 0]
eD2_B1 = [x for (x, mask1, mask2) in zip(eD2_B, eD1_B, aMask_B) if mask1 != 0 and mask2 == 0]

# Exponential

ratio_e_A = map(truediv, eE2_A1, eE1_A1)
angle_e_A = [0.5*orientation(x) for x in ratio_e_A]

ratio_e_B = map(truediv, eE2_B1, eE1_B1)
angle_e_B = [0.5*orientation(x) for x in ratio_e_B]

plt.hist(angle_e_A, bins=60, color='r')
plt.hist(angle_e_B, bins=60, color='b', alpha=0.7)
plt.xlabel("Angle (degrees)")
plt.title("Angle of Orientation of Ellipticity Histogram (Exponential)")
plt.show()

# de Vaucouleurs
ratio_d_A = map(truediv, eD2_A1, eD1_A1)
angle_d_A = [0.5*orientation(x) for x in ratio_d_A]

ratio_d_B = map(truediv, eD2_B1, eD1_B1)
angle_d_B = [0.5*orientation(x) for x in ratio_d_B]

plt.hist(angle_d_A, bins=60, color='r')
plt.hist(angle_d_B, bins=60, color='b', alpha=0.7)
plt.xlabel("Angle (degrees)")
plt.title("Angle of Orientation of Ellipticity Histogram (de Vaucouleurs)")
plt.show()


# Mean, Median, and Moment (Skewness)
a = numpy.mean(angle_e_A)
b = numpy.mean(angle_e_B)
c = numpy.mean(angle_d_A)
d = numpy.mean(angle_d_B)

print('Means (EXP A/B, DEV A/B):', a, b, c, d)

x = moment(angle_e_A, moment=3)
y = moment(angle_e_B, moment=3)
z = moment(angle_d_A, moment=3)
w = moment(angle_d_B, moment=3)

print('Moment (EXP A/B, DEV A/B):', x, y, z, w)
