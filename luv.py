import numpy as np
from astropy.io import fits
import math
import matplotlib.pyplot as plt
from cosmolopy import fidcosmo, magnitudes, cc, cd
import cosmolopy.distance as cd
import scipy
from scipy.integrate import simps
c = 299792.458
cosmo = {'omega_M_0': 0.3, 'omega_lambda_0': 0.7, 'omega_k_0': 0.0, 'h': 0.70}
z=2.8
dl = cd.luminosity_distance(z, **cosmo) #Mpc
da=dl/((1.+z)**2)*1e3 #kpc

img, ih = fits.getdata('img_c_4.fits', header=True)  # f_lamda: erg s-1 cm-2 A-1

# L3
center_x=23
center_y=27
w_x = 1 #aperture
w_y = 2

def cut(wave_l,wave_u,image):
    l=int((wave_l-4075)/0.5)
    u=int((wave_u-4075)/0.5)
    band = np.sum(image[l:u,:,:], axis=0)*0.5
    return band

band = cut(4610,4640,img)-cut(4550,4580,img)  # erg s-1 cm-2; flux integration
flag = np.zeros(band.shape)  # flag for aperture
flag[(center_y-w_y)+1:(center_y+w_y)+1,(center_x-w_x)+1:(center_x+w_x)+1] = 1
L_lya_L3= 4*3.14*(dl*3.086*10**24)**2*np.sum(band[flag==1])*(1.+z)**3  # luminousity lamda lya: erg s-1;  aperture integration

s=da**2*abs(len(flag[flag==1])*ih['CD2_2']*ih['CD1_1'])*(2*3.14/360)**2 # area of L3: kpc^2
R=7.5/3600*2*3.14/360*da   # ditance between L1 and L3: kpc
d_sigma=s/(R**2) # solid angle of L3
frac = d_sigma/(4*3.14)
L_uv=L_lya_L3/0.6/frac
print L_uv




