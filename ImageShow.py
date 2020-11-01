import astropy.io.fits as fits
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


def openImage(file,inverted):
    fitsStars = fits.open(file)
    imgStars = fitsStars[0].data.astype(np.float64)
    showImage(imgStars, "Image", inverted)


openImage("Images/FitsImages/L_2019-04-08_22-59-51_c.fits",True)
