import matplotlib.pyplot as plt
from matplotlib import cm
import astropy as ap
from astropy.io import fits
import numpy as np
import argparse
import os

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--image", required = True)
ap.add_argument("-c", "--corr", required = True)
args = vars(ap.parse_args())

fitsImg = fits.open(args['image'])
img = fitsImg[0].data
mx = max([max(img[i]) for i in range(0,len(img))])/18
mn = min([min(img[i]) for i in range(0,len(img))])
plt.imshow(img, cmap=cm.gray, vmin = mn, vmax = mx)
plt.show()


fitsCorr = fits.open(args['corr'])
corr = fitsCorr[1].data

xs = []
ys = []
k=0
for item in corr:
    x = int(round(item[4]))
    y = int(round(item[5]))
    #sl = img[x-20:x+20,y-20:y+20]
    #if k==1:
    #    print(sl)
    #result = np.where(sl == np.amax(sl))
    #x+=result[0]
    #y+=result[1]
    for i in range(-10,10):
        for j in range(-10,10):
            try:
                img[y+i][x+j]=0
            except:
                pass
    #k+=1
    
plt.imshow(img, cmap=cm.gray, vmin = mn, vmax = mx)
plt.show()
