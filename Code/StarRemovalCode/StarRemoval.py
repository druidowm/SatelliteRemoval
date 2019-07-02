import matplotlib.pyplot as plt
import pyplate
import astropy as ap
from astropy.io import fits
import numpy as np
#import cv2
import argparse
import os

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required = True)
ap.add_argument("-o", "--output", required = True)
args = vars(ap.parse_args())

psPipe = pyplate.pipeline.PlatePipeline()
psPipe.assign_conf("PlateSolveConfig.txt")

psPipe.single_image(args["input"])

print("here")
fitsImg = fits.open(os.path.join("../../Images/FitsImages",args["input"]))#Readonly
print("here")
fitsHeader = fitsImg.info()
print(fitsHeader)
print(fitsImg[0].data)
