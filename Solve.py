import os
import TrailRemovalCode.RemoveTrail as RT
import StarRemovalCode.StarRemoval as SR
import TrailRemovalCode.GenerateTrail as GT
import TrailRemovalCode.RemoveBackground as RB
from photutils.background import MMMBackground
from TrailRemovalCode.PiecewiseCubic import PiecewiseCubic as PC
import numpy as np
import astropy as ap
from astropy.io import fits
from matplotlib import cm
from matplotlib import pyplot as plt
from photutils import Background2D

rescale = 200
rescale2 = 30

def solveWithClient(imgPath, outPath, trailSize, showProgress = True, showNoStars = False, showTrail = False, showNoTrail = False, showNoTrailWithTrail = False, inverted = False):
    if showProgress:
        print("\n\n\n\n")

    fitsStars = fits.open(imgPath)
    imgStars = fitsStars[0].data.astype(np.float64)

    imgNoBackground = imgStars - Background2D(imgStars,100).background

    SR.showImage(imgNoBackground, "Image with Background Subtracted", inverted)

    SR.showImage(imgStars, "Original Image", inverted)

    imgNoStars,gX,gY,gR,starX,starY,starR,mags,hotX,hotY = SR.findAndRemoveStars(np.copy(imgNoBackground), imgPath, showProgress = showProgress, inverted = inverted)


    imgNoTrail = RT.findRemoveTrail(imgNoStars, imgNoBackground, imgStars, starX, starY, starR, gX, gY, gR, hotX, hotY, trailSize, showProgress, showNoStars, showTrail, showNoTrail, showNoTrailWithTrail, inverted)
    

    fits.writeto(outPath, imgNoTrail, header = fitsStars[0].header)

    if showProgress:
        print("Satellites B gone!")

def batchSolveWithClient(directory, outDir, trailSize, showProgress = True, showNoStars = False, showTrail = False, showNoTrail = False, showNoTrailWithTrail = False, inverted = False):
    images = os.scandir(directory)

    i = 1
    for image in images:
        if image.name != ".DS_Store" and not image.is_dir():
            print("\n\n\n\n\n\n\n\n")
            print("Removing trail for image " + str(i)+" "+str(directory)+"/"+str(image.name))
            dotPos = (image.name).rfind('.')
            outPath = outDir+"/"+image.name[0:dotPos]+"_no_trail"+image.name[dotPos:]
            solveWithClient(directory+"/"+image.name, outPath, trailSize, showProgress, showNoStars, showTrail, showNoTrail, showNoTrailWithTrail, inverted)
            i+= 1


def solveTest(imgPath, outPath, trailSize, p1, p2, showTrail = False, showNoTrail = False, showNoTrailWithTrail = False, inverted = False):

    fitsStars = fits.open(imgPath)
    imgStars = fitsStars[0].data.astype(np.float64)

    imgNoBackground = imgStars - Background2D(imgStars,100).background

    imgNoStars,gX,gY,gR,starX,starY,starR,mags,hotX,hotY  = SR.findAndRemoveStars(np.copy(imgNoBackground), imgPath, showProgress = True, inverted = inverted)
    
    imgNoStarsPSF, starX, starY, starR = SR.removeStarsPSF(np.copy(imgNoBackground), starX, starY, starR, mags, gX, gY, gR, p1, p2, 2*np.mean(starR), showProgress = True, inverted = inverted)

    fit= PC(imgNoBackground)

    imgNoTrail = RT.removeTrail(imgNoBackground, imgStars, True, p1, p2, starX,starY,starR,hotX,hotY, fit, trailSize)
    
    RT.showImage(imgNoTrail, "Image Without Trail", inverted)
    
    fits.writeto(outPath, imgNoTrail, header = fitsStars[0].header)

    
def testTrailFinder(xScale, yScale):
    for i in range(0,10):
        trail = GT.generateTrailVertical(np.random.random()*xScale, np.random.random(), np.random.random(), np.random.random(), 50*np.random.random(), xScale, yScale)
        imgNoTrail = RT.findRemoveTrail(trail, trail, 200, False, True, True, True, False, False)
