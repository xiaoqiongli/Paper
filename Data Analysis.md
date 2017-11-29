# Data Analysis




## Sky subtract
To import our functions

	import sub_sky as ss
	
To create a a new folder 'analysis'

	ss.subsky()
	ls #check
	
Subtract sky:

Be careful to check *_mask.fits whether mask all the pixels of source or not.

	ss.sub_median('kb171020_00041_icubes.fits', 'kb171020_00040_icubes.fits', '40_41', 2, cut_ch=300)
	ss.sub_median('kb171020_00041_icubes.fits', 'kb171020_00042_icubes.fits', '42_41', 2, cut_ch=300)
	ss.sub_median('kb171020_00044_icubes.fits', 'kb171020_00043_icubes.fits', '43_44', 2, cut_ch=300)
	ss.sub_median('kb171020_00046_icubes.fits', 'kb171020_00045_icubes.fits', '45_46', 2, cut_ch=300)

Subtract continuum source, not good, it's optional. 


	ss.psf('40_41', dpx=8, dpy=10, scale='integrated')

Anyway, we can subtract continuum when we plot the image (i.e. line emission subtract the continue: Lya img[4610A:4640A]-img[4550A:4580A])
	
##Combine

Combine images as the peak position

	ss.combine('40_41.fits','42_41.fits','img_c_2.fits')
	ss.combine('img_c_2.fits','43_44.fits','img_c_3.fits')
	ss.combine('img_c_3.fits','45_46.fits','img_c_4.fits')


