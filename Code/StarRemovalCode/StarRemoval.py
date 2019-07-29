import matplotlib.pyplot as plt
from matplotlib import cm
import astropy as ap
from astropy.io import fits
from astropy import wcs
import astropy.coordinates as coord
import astropy.units as u
from astroquery.vizier import Vizier
#from astroquery.gaia import Gaia
import numpy as np
import math
import argparse
import os

ap = argparse.ArgumentParser()
#The image that needs to be processed
ap.add_argument("-i", "--image", required = True)
ap.add_argument("-c", "--corr", required = True)
args = vars(ap.parse_args())

fitsImg = fits.open(args['image'])
img = fitsImg[0].data
mx = max([max(img[i]) for i in range(0,len(img))])/18 #8 for Ryan's image #18 for test.fits
mn = min([min(img[i]) for i in range(0,len(img))])
plt.imshow(img, cmap=cm.gray, vmin = mn, vmax = mx)
plt.show()

#w = wcs.WCS(args["image"])
#topLeft = w.all_pix2world(0, 0, 0)
#botRight = w.all_pix2world(len(img), len(img[0]), 0)
#center = [round((topLeft[0]+botRight[0])/2),round((topLeft[1]+botRight[1])/2)]#

#v = Vizier(columns=['_RAJ2000', '_DEJ2000', '_VTmag'],
#           column_filters={"VTmag":">15"})

#v.ROW_LIMIT = -1


#objects = v.query_region(coord.SkyCoord(ra = center[0], dec = center[1], unit = (u.deg,u.deg), frame = 'icrs'),
#                              radius = str(math.sqrt((topLeft[0]-botRight[0])**2+(topLeft[1]-botRight[1])**2))+"d")

#gaia = objects["I/337/gaia"]
#print(gaia.columns)
#print(gaia)
#print(gaia[0]["_RAJ2000"])

#objects = Gaia.query_object_async(coordinate=coord.SkyCoord(ra = center[0], dec = center[1], unit = (u.deg,u.deg), frame = 'icrs'), width=str(topLeft[0]-botRight[0])+"d", height=str(topLeft[1]-botRight[1])+"d")

#print(objects)
mean = np.mean(img)
#print(mean)

fitsCorr = fits.open(args['corr'])
corr = fitsCorr[1].data

xs = []
ys = []
for item in corr:
    x = int(round(item[0]))
    y = int(round(item[1]))

    inStar = True
    r=1

    while inStar:
        for theta in range(0,8,1):
            boundX = round(x+r*math.cos(theta*math.pi/4))
            boundY = round(y+r*math.sin(theta*math.pi/4))


            try:
                if img[boundY][boundX] <= mean:
                    inStar=False
                    break
            except:
                pass

        r+=1

    r-=1
    for i in range(-r,r):
        h = round(math.sqrt(r**2-i**2))
        for j in range(-h,h):
            try:
                img[y+i][x+j]=mean
            except:
                  pass
    
plt.imshow(img, cmap=cm.gray, vmin = mn, vmax = mx)
plt.show()
