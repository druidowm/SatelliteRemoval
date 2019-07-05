import matplotlib.pyplot as plt
import astropy as ap
from astropy.io import fits
import numpy as np
import argparse
import os

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required = True)
args = vars(ap.parse_args())

fitsImg = fits.open(args['input'])#Readonly
fitsHeader = fitsImg.info()
print(fitsHeader)
for item in fitsImg:
    print(item.data)
