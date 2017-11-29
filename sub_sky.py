""" Python Vesion: 2.7.13
    Functions to subtract sky:
"""

def subsky():
    import os
    newpath = r'analysis'
    if not os.path.exists(newpath):
        os.makedirs(newpath)


def sub_median(sky_file, img_file, output, mask_limit, cut_ch=300):
    """ This function is used to subtract sky median

    Parameters
    ----------
    sky_file
    img_file
    mask_limit : mask limit
    cut_ch     : cut bad channel number

    Eaxmple:
    sub_median('kb171021_00083_icuber.fits', 'kb171021_00082_icuber.fits', 'img', 40, cut_ch=300)

    Returns
    -------
    mask.fits  : masked image (mask the bad pixels and qso)
    img.fits   : sky median subtracted image
    """
    from astropy.io import fits
    import numpy as np
    from astropy.stats import sigma_clip
    sky, sh = fits.getdata(sky_file, header=True)
    img, ih = fits.getdata(img_file, header=True)

    ### begin to mask qso on sky
    # define final mask image
    mask_sum = np.zeros(sky[1,:,:].shape)

    ## looping over the layers
    for j in np.arange(cut_ch, (len(img[:,1,1])-cut_ch), 1):
        # iterated 10 times, rejecting points by > 3 sigma
        filtered_data = sigma_clip(sky[j,:,:], sigma=3, iters=10)
        sky_value = np.median(filtered_data.data[filtered_data.mask == False])
        # subtract the sky median value
        sky[j,:,:]-=sky_value
        img[j,:,:]-=sky_value
        # record the masked counts
        mask = np.zeros(sky[1,:,:].shape)
        mask[filtered_data.mask == True] = 1
        mask_sum+=mask

    # use mask_limit to avoid masking the good pixels (noise)
    mask_sum[mask_sum<mask_limit]=0
    mask_sum[mask_sum>=mask_limit]=1

    fits.writeto('analysis/'+output+'_mask.fits',mask_sum, overwrite=True)
    fits.writeto('analysis/'+output+'.fits',img, ih, overwrite=True)



def sub_psf(img_file, dpx, dpy, scale):
    """ This function is used to subtract PSF for bright point source

    Parameters
    ----------
    img_file  :   default 'img.fits'
    dpx       :   cut PSF region in x
    dpy       :   cut PSF region in y
    scale     :   sclae method: 'peak', 'p', 'integrated', 'i'

    Eaxmple:
    psf('img.fits', dpx=8, dpy=10, scale='integrated')

    Returns
    -------
    img_psf.fits  : psf subtracted image
    """
    from astropy.io import fits
    import numpy as np
    img, ih = fits.getdata('analysis/'+img_file+'.fits', header=True)

    ### PSF reduced
    img_median= np.median(img, axis=0)

    # find peak location
    loc = np.where(img_median > 0.9*np.amax(img_median))
    px= int(round(np.mean(loc[1])))
    py= int(round(np.mean(loc[0])))

    a = np.zeros(img[:,1,1].shape)
    psf = np.zeros(img.shape)

    ## scale to each channel
    if scale not in ['peak', 'p', 'integrated', 'i']:
        raise IOError('Need input the scale method! scale = peak / integrated')
    # peak value to scale
    if scale in ['peak', 'p']:
        peak_value = np.max(img_median)
        a=img[:,py,px]/peak_value

    # cut the region
    for j in np.arange(0, len(img[:,1,1]), 1):
        # integrated to scale
        if scale in ['integrated', 'i']:
            integrated_value = np.sum(img[j,(py-dpy):(py+dpy),(px-dpx):(px+dpx)])
            a[j]=integrated_value/np.sum(img_median[(py-dpy):(py+dpy),(px-dpx):(px+dpx)])
        psf[j,(py-dpy):(py+dpy),(px-dpx):(px+dpx)]=img_median[(py-dpy):(py+dpy),(px-dpx):(px+dpx)]*a[j]

    # reduce psf
    img-=psf

    fits.writeto('analysis/'+img_file+'_psf.fits',img, ih)


def cube_shift(image,shift_x,shift_y):
    '''
      This function is used to offset a datacube
     INPUT:     - img_file
                - offset x, offset y (in pixels)
                - offset can be determined using find_center.py
     OUTPUT:    - shifted datacube for combining
     INPUT EXAMPLE:

       img, ih = fits.getdata("kb171021_00077_icuber.fits", header=True)
       img_shift= cube_shift(img,shift_x,shift_y)
     EXAMPLE 2:

       img1, ih1 = fits.getdata("kb171021_00077_icuber.fits", header=True)
       img2, ih2 = fits.getdata("kb171021_00082_icuber.fits", header=True)
       center1=find_center(img1)
       center2=find_center(img2)
       shift_x= center2[0]-center1[0]
       shift_y= center2[1]-center1[1]
       img2_shift= cube_shift(img2,shift_x,shift_y)
       center_new= find_center(img2_shift)
    Then One can see new center for img2 is the same as that of img1. Ready for combining
     MODIFICATION HISTORY:
       2017-11-02 Initial version (ZC)
    '''
    import numpy as np
    from scipy import interpolate
    image_sh= np.copy(image)  # construct a new datacube called image_sh (image_shift)

    for i in np.arange(0, len(image[:,1,1]), 1):  # loop layers
        # offset each layer using shift_x, shift_y
        image_sh[i,:,:]= np.roll(image[i,:,:], -int(round(shift_x)), axis=0)
        image_sh[i,:,:]= np.roll(image_sh[i,:,:], -int(round(shift_y)), axis=1)

    return image_sh

def find_center(image):
    '''
    Find the center of QSOs for each target,
    we will return the center of the QSO;
    and later, we will use the center to
    calculate the offset.  (Z. Cai)
    This function is used to determine the center of the QSO for each datacube
    INPUT EXAMPLE:

    img= fits.getdata("kb171021_00092_icuber.fits", header=True)
    centers= find_center(img)
    MODIFICATION HISTORY:
    2017-11-01 Initial version (ZC)
    '''
    import numpy as np
    from scipy import interpolate
    # construct a 2D broadband image using median stacking along the wavelength
    cut_ch=300
    img_c=image[cut_ch:(len(image[:,1,1])-cut_ch),:,:]
    sci_med=np.median(img_c, axis=0)

    # interpolate to a finer grid:
    x_ori= np.linspace(0, len(sci_med[:,0])-1, len(sci_med[:,0])) # 0--23 (24 pixels)
    y_ori= np.linspace(0, len(sci_med[0,:])-1, len(sci_med[0,:])) # 0--69 (70 pixels)
    x_new= np.linspace(0, len(sci_med[:,0])-1, 3.*len(sci_med[:,0])-2)  # new grid, 3x pixel number
    y_new= np.linspace(0, len(sci_med[0,:])-1, 3.*len(sci_med[0,:])-2)  # new grid, 3x pixel number

    f = interpolate.interp2d(y_ori, x_ori, sci_med,kind='cubic')
    sci_interp = f(y_new,x_new)

    # find center using pixels > 0.93* the peak value:
    cut_x=30
    cut_y=4
    a=sci_med[cut_x:len(sci_med[:,0])-cut_x,cut_y:len(sci_med[0,:])-cut_y]
    indices = np.where(a > 0.93*np.amax(a))
    #indices_interp = np.where(sci_interp > 0.93*np.amax(sci_interp))

    # center using the original pixel grid
    x_mean0= int(np.mean(indices[0])+cut_x)
    y_mean0= int(np.mean(indices[1])+cut_y)

    # center using the new finer pixel grid (1/3 x original pixel scale)
    #x_mean_interp= np.mean(indices_interp[0])*1/3.
    #y_mean_interp= np.mean(indices_interp[1])*1/3.

    # return center determined using both the original and new pixel grid.
    #center= [x_mean0, y_mean0, x_mean_interp,y_mean_interp]
    center= [x_mean0, y_mean0]
    #fits.writeto('peak.fits',sci_med, ih1)
    return center

def combine(image1,image2,output):
    '''
    '''
    from astropy.io import fits
    import numpy as np
    img1, ih1 = fits.getdata('analysis/'+image1, header=True)
    img2, ih2 = fits.getdata('analysis/'+image2, header=True)

    center1=find_center(img1)
    center2=find_center(img2)

    shift_x= center2[0]-center1[0]
    shift_y= center2[1]-center1[1]

    img2_shift= cube_shift(img2,shift_x,shift_y)

    img_comb=img1+img2_shift
    center=find_center(img_comb)
    fits.writeto('analysis/'+output,img_comb, ih1, overwrite=True)