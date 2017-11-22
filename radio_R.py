import numpy as np
import cosmolopy.distance as cd

cosmo = {'omega_M_0': 0.3, 'omega_lambda_0': 0.7, 'omega_k_0': 0.0, 'h': 0.70}

z=2.803
m_B=23.88

dl = cd.luminosity_distance(z, **cosmo) #Mpc
m_B=m_B-0.142 # correct extinction from D. Schlegel et al.(1998ApJ...500..525S)
M_B=m_B-5*(np.log10(dl*10**6)-1)

f_B=4.26*10**(-20-0.4*m_B)   #https://en.wikipedia.org/wiki/Apparent_magnitude
f_B2=10**(-0.4*(m_B+48.574)) #https://en.wikipedia.org/wiki/Jansky
f_r=(5**(-1.2))/(1.4**(-1.2))*526*10**(-23-6) # erg s-1 cm-2 Hz-1   ##### 0<a<1 why 1.2 ? ######
f_r2=(5**(-1.2))/(8.7**(-1.2))*57*10**(-23-6) # upper limit


print f_r/f_B



