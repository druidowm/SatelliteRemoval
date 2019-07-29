import matplotlib.pyplot as plt
import astropy as ap
from astropy.io import fits
import numpy as np
import argparse
import os

def printData(d):
    if type(d) is fits.ImageHDU:
        print(1)
        print(d.data)
    elif type(d) is fits.HDUList:
        print(2)
        for item in d:
            printData(item)
    elif type(d) is fits.PrimaryHDU:
        print(3)
        print(d.data)

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required = True)
args = vars(ap.parse_args())

fitsImg = fits.open(args['input'])#Readonly
fitsHeader = fitsImg.info()
print(fitsHeader)
printData(fitsImg)

        
