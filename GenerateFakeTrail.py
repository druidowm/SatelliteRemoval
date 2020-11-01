import TrailRemovalCode.GenerateTrail as GT
import astropy as ap
from astropy.io import fits
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np

rescale = 200
rescale2 = 30

def showImage(img, title, inverted):
    mx = img.max()
    mn = img.min()
    median = np.median(img)

    if inverted:
        plt.imshow(-img, cmap=cm.gray, vmin = -(mx+rescale*median)/(rescale+1), vmax = -(mn+rescale2*median)/(rescale2+1))
    else:
        plt.imshow(img, cmap=cm.gray, vmin = (mn+rescale2*median)/(rescale2+1), vmax = (mx+rescale*median)/(rescale+1))

    plt.title(title)
    
    plt.show()

inFile = "Images/FitsImages/2019-12-22_04-23-01__-30.00_60.00s_185_c.fits"
outFile = "Images/FitsImages/2019-12-22_04-23-01__-30.00_60.00s_185_c_double_trail.fits"

fitsImg = fits.open(inFile)
img = fitsImg[0].data.astype(np.float64)
for i in range(0,10):
    trail = GT.generateTrailHorizontal(4*np.random.random()-2, 3000*np.random.random(), 0.0003*np.random.random()-0.000005, 0.03*np.random.random()-0.0005, 50*np.random.random()+50, 5*np.random.random(), img.shape[1], img.shape[0])

    imgTrail = img + trail

    showImage(imgTrail,"Image with Artificial Trail", inverted = False)


    fits.writeto(outFile+str(i+1), imgTrail, header = fitsImg[0].header)
