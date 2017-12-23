import numpy as np
import math
import matplotlib.pyplot as plt
from cosmolopy import fidcosmo, magnitudes, cc, cd
import cosmolopy.distance as cd
import scipy
from scipy.integrate import simps
c = 2.99792458e18 #Angs/s
cosmo = {'omega_M_0': 0.3, 'omega_lambda_0': 0.7, 'omega_k_0': 0.0, 'h': 0.70}
z=2.8

dl = cd.luminosity_distance(z, **cosmo) #Mpc

# Zero Magnitude Flux Density:
wiseband=1     # WISE band
wiseband=wiseband-1
w_lamd=np.array([   3.4,     4.6,      12,     22])*1e4
m_w   =np.array([16.478,  15.267,  12.327,  9.311])
m_werr=np.array([ 0.066,   0.085,   0.342,  0.491])
def funfv(m,wiseband): # Eq 1 in http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#example
    fv0=np.array([309.540,171.787,31.674,8.363])
    fv0=fv0[wiseband]
    fv=fv0*10**(-m/2.5)
    return fv
fv=funfv(m_w[wiseband],wiseband)*1e-23  # erg s-1 cm-2 Hz-1
fverr=abs(funfv(m_w[wiseband],wiseband)-funfv(m_w[wiseband]-m_werr[wiseband],wiseband))*1e-23  # erg s-1 cm-2 Hz-1
vfv=fv*c/w_lamd[wiseband]


# Type II QSO template: http://www.iasf-milano.inaf.it/~polletta/templates/swire_templates.html
data=[]
r = open('type-II.sed','r')
a= r.readlines()
for i in a:
    b = i.strip().split()
    data.append(b)
data=np.array(data,dtype=float)
r.close()
l = data[:,0]
norm_fl = data[:,1]  # norm f_lamda

fl_w = vfv/w_lamd[wiseband]*(1.+z)**3 # rest frame F_lamda at WISE
fl = norm_fl*fl_w/norm_fl[abs(w_lamd[wiseband]/(1.+z)-l)==np.min(abs(w_lamd[wiseband]/(1.+z)-l))] #erg-1 s-1 cm-2
Ll = 4*3.14*(dl*3.086*10**24)**2*fl # luminousity lamda
L_iso = scipy.integrate.trapz(Ll,l) #erg-1 s-1; intergration for all wave



### NOT GOOD METHOD:
###infared bolometric correction   (https://arxiv.org/pdf/1207.2124.pdf)
###WISE1: 3.4um 1.5um
###WISE2: 4.6um 1.5um
###WISE3: 12um  3um
###WISE4: 22um  7um
##l_llamd= 4*3.14*(dl*3.086*10**24)**2*vfv*(1.+z)**2  # rest frame F_lamda at WISE
##def fun_liso(l_llamd,wiseband):
##    k=np.array([0.82,  0.82,  0.92,  0.88])
##    b=np.array([8.98,  8.98,  4.54,  1.85])
##    l_iso=10**(b[wiseband] + k[wiseband] * np.log10(l_llamd))
##    return l_iso
##
##L_iso=fun_liso(l_llamd,wiseband)

# additional correction to bolometric luminosity for the assumption of isotropy
L = 0.75 * L_iso

print L



