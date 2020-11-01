import astropy
import photutils
from photutils.psf import DAOPhotPSFPhotometry as DAOP
from photutils.psf import IntegratedGaussianPRF as PRF
from photutils.background import MMMBackground
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import astropy.io.fits as fits
from matplotlib.colors import LogNorm
from photutils.datasets import (make_random_gaussians_table,
                                make_noise_image,
                                make_gaussian_sources_image)
from photutils.psf import (IterativelySubtractedPSFPhotometry,
                           BasicPSFPhotometry)
from photutils import MMMBackground
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.detection import DAOStarFinder
from photutils.detection import IRAFStarFinder
from photutils.utils import calc_total_error
from astropy.table import Table
from astropy.modeling.fitting import LevMarLSQFitter


#x = np.arange(0,50)
#y = np.arange(0,50)
#x,y = np.meshgrid(x,y)
#image = x*0.0
#for i in range(0,20):
#    z = -((x-50*np.random.random())**2+(y-50*np.random.random())**2)/7
#    image += 50*np.random.random()*np.exp(z)

#plt.imshow(image)
#plt.show()

fitsStars = fits.open("Images/FitsImages/L_2019-04-08_22-59-51_c.fits")
imgStars = fitsStars[0].data.astype(np.float64)



#bkg = MMMBackground()
#background = bkg(imgStars)
#gaussian_prf = PRF()
#gaussian_prf.sigma.fixed = False
#photTester = DAOP(8,background,5,gaussian_prf,(11,11))

daogroup = DAOGroup(crit_separation=8)
mmm_bkg = MMMBackground()
iraffind = DAOStarFinder(threshold=2*mmm_bkg(imgStars), fwhm=4.5)
fitter = LevMarLSQFitter()
gaussian_prf = IntegratedGaussianPRF(sigma=2.05)
gaussian_prf.sigma.fixed = False
photTester = IterativelySubtractedPSFPhotometry(finder=iraffind,
                                                   group_maker=daogroup,
                                                   bkg_estimator=mmm_bkg,
                                                   psf_model=gaussian_prf,
                                                   fitter=fitter,
                                                   fitshape=(11, 11),
                                                   niters=2)

photResults = photTester(imgStars)
print(photResults['x_fit','y_fit','flux_fit'])

finalImg = photTester.get_residual_image()

rescale = 200
rescale2 = 30

mx1 = finalImg.max()
mn1 = finalImg.min()
median1 = np.median(finalImg)

mx2 = imgStars.max()
mn2 = imgStars.min()
median2 = np.median(imgStars)

print(median1)
print((mx1+median1)/(2))
print((mn1+median1)/(2))

fig,axes = plt.subplots(1,2,sharex=True,sharey=True)
axes[0].imshow(-finalImg, cmap=cm.gray)#, vmin = -(mx1+100000*rescale*median1)/(100000*rescale+1), vmax = -(mn1+100000*rescale2*median1)/(100000*rescale2+1))
#axes[0].plot(photResults['x_fit'],photResults['y_fit'],"ro")
axes[1].imshow(-imgStars, cmap=cm.gray)#, vmin = -(mx2+rescale*median2)/(rescale+1), vmax = -(mn2+rescale2*median2)/(rescale2+1))
plt.show()

