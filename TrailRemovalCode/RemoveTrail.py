import TrailRemovalCode.DetectSatellite as DS
import TrailRemovalCode.GenerateTrail as GT
import TrailRemovalCode.TrailFit as TF
from TrailRemovalCode.PiecewiseCubic import PiecewiseCubic as PC
from photutils.background import MMMBackground

from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import Polynomial as P
import math
from scipy.optimize import curve_fit as cf

rescale = 200
rescale2 = 30

def showImageWithTrail(img, title, p1, p2, inverted, offset = None):
    mx = img.max()
    mn = img.min()
    median = np.median(img)

    if inverted:
        plt.imshow(-img, cmap=cm.gray, vmin = -(mx+rescale*median)/(rescale+1), vmax = -(mn+rescale2*median)/(rescale2+1))
    else:
        plt.imshow(img, cmap=cm.gray, vmin = (mn+rescale2*median)/(rescale2+1), vmax = (mx+rescale*median)/(rescale+1))

    plt.plot((p1[0],p2[0]),(p1[1],p2[1]))
    if offset != None:
        m,b = TF.getLineParams(p1, p2)
        pars = TF.getParLine(img, p1, p2, m, b)
        perps = offset(pars)
        x,y= GT.getXY(pars,perps,m,b)
        plt.plot(x,y)
    plt.show()

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

def getLineParams(p1, p2):
    m = (p2[1]-p1[1])/(p2[0]-p1[0])
    b = p1[1]-m*p1[0]
    return (m,b)

def removeTrail(imgNoBackground, img, showProgress, p1, p2, starX, starY, starR, gX, gY, gR, hotX, hotY, imgPC, trailSize):
    if abs(p1[0]-p2[0]) < 0.01:
        return removeTrailVertical(imgNoBackground, img, showProgress, p1[0], trailSize)
    else:
        m,b = getLineParams(p1, p2)
        return removeTrailHorizontal(imgNoBackground, img, showProgress, m, b, p1, p2, starX, starY, starR, gX, gY, gR, hotX, hotY, trailSize, imgPC)

def removeTrailHorizontal(imgNoBackground, img, showProgress, m, b, p1, p2, starX, starY, starR, gX, gY, gR, hotX, hotY, trailSize, PC):
    if showProgress:
        print("Fitting curve brightness...")

    sDev, brightness, offset = TF.findFitHorizontal(imgNoBackground, p1, p2, starX, starY, starR, gX, gY, gR, hotX, hotY, trailSize, PC)
    
    if showProgress:
        print("Brightness fitted\n")
        print("Removing trail...")
    trail = GT.generateTrailHorizontal(m, b, brightness, sDev, img.shape[1], img.shape[0], offset)
    if showProgress:
        print("Trail removed!")
    imgNoTrail = img - trail
    showImageWithTrail(imgNoTrail, "Image Without Trail", p1, p2, False, offset)
    mx = img.max()
    mn = img.min()
    median = np.median(img)

    fig,axes = plt.subplots(1,2,True,True)

    axes[0].imshow(-img, cmap=cm.gray, vmin = -(mx+rescale*median)/(rescale+1), vmax = -(mn+rescale2*median)/(rescale2+1))
    axes[1].imshow(-imgNoTrail, cmap=cm.gray, vmin = -(mx+rescale*median)/(rescale+1), vmax = -(mn+rescale2*median)/(rescale2+1))
    plt.show()
    
    return imgNoTrail

def removeTrailVertical(imgNoBackground, img, showProgress, x, trailSize):
    if showProgress:
        print("Fitting curve brightness...")
    b1, b2, b3, sDev = TF.findFitVertical(imgNoBackground, x, trailSize)
    if showProgress:
        print("Brightness fitted\n")
        print("Removing trail...")
    print(sDev)
    trail = GT.generateTrailVertical(x, b1, b2, b3, sDev, img.shape[1], img.shape[0])
    if showProgress:
        print("Trail removed!")
    return img - trail

def findRemoveTrail(imgNoStars, imgNoBackground, imgStars, starX, starY, starR, gX, gY, gR, hotX, hotY, trailSize, showProgress, showNoStars = False, showTrail = False, showNoTrail = False, showNoTrailWithTrail = False, inverted = False):
    if showNoStars:
        showImage(imgNoStars, "Image Without Stars", inverted)
    if showProgress:
        print("Searching for satellite trail...")
    p1,p2 = DS.findSatelliteScale(imgNoStars, debug = True, imgStars = imgNoBackground)

    imgPC = PC(imgNoStars)

    if showProgress:
        print("Satellite trail found\n")
        
    if showTrail:
        showImageWithTrail(imgNoStars, "Image With Detected Trail", p1, p2, inverted)

    imgNoTrail = removeTrail(imgNoBackground, imgStars, showProgress, p1, p2, starX, starY, starR, gX, gY, gR, hotX, hotY, imgPC, trailSize)

    if showNoTrail:
        if showNoTrailWithTrail:
            showImageWithTrail(imgNoTrail, "Image Without Trail", p1, p2, inverted)

        else:
            showImage(imgNoTrail, "Image Without Trail", inverted)
    
    return imgNoTrail
