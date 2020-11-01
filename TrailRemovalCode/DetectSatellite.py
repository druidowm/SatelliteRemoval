import math
import numpy as np
import scipy
import scipy.stats
from scipy.optimize import curve_fit as cf
import matplotlib.pyplot as plt
import cv2
import TrailRemovalCode.RemoveTrail as RT
from TrailRemovalCode.PiecewiseCubic import PiecewiseCubic as PC


def frange(start, stop, step):
    x = start
    while x<=stop:
        yield x
        x += step

def findSatelliteScale(imdata, fit=None, debug = False, iter = 1, imgStars = None):
    if fit == None:
        #if type(imgStars) == type(None):
        fit = PC(imdata) 
        #else:
        #    fit = PC(imgStars)

    if imdata.shape[0]*imdata.shape[1] <= 90000:
        return findSatellite(imdata, fit)
    
    newX = imdata.shape[1]//4
    newY = imdata.shape[0]//4
    xScale = imdata.shape[1]/newX
    yScale = imdata.shape[0]/newY

    imgScaledDown = cv2.resize(imdata, (newX,newY), cv2.INTER_AREA)
    if debug:
        RT.showImage(imgScaledDown, "Image With Trail", False)

    fit2 = PC(imgScaledDown)
    p1,p2,filler = findSatelliteScale(imgScaledDown, fit2, debug, iter+1)

    if debug:
        RT.showImageWithTrail(imgScaledDown, "Image With Trail", (p1[0],p1[1]), (p2[0],p2[1]), False)

    step = 1
    p1Tests = findTestLines(p1, xScale, yScale, imdata.shape[1], imdata.shape[0], newX, newY, step)
    p2Tests = findTestLines(p2, xScale, yScale, imdata.shape[1], imdata.shape[0], newX, newY, step)

    maxCount = -1
    maxP1 = (0.0,0.0)
    maxP2 = (0.0,0.0)

    minLen = max(imdata.shape[0],imdata.shape[1])/4

    for i in range(0,len(p1Tests)):
        for j in range(0,len(p2Tests)):
            count = searchLine(fit, p1Tests[i], p2Tests[j], minLen)
            if count > maxCount:
                maxCount = count
                maxP1 = p1Tests[i]
                maxP2 = p2Tests[j]

    return (maxP1,maxP2,fit)


def findTestLines(p, xScale, yScale, oldX, oldY, newX, newY, step):
    tests = []
    if p[0] == 0:
        if p[1] == 0:
            for i in frange(0,math.ceil(xScale),step):
                tests.append((i,0))
            for i in frange(1,math.ceil(yScale),step):
                tests.append((0,i))
        
        elif p[1] == newY - 1:
            for i in frange(0,math.ceil(xScale),step):
                tests.append((i,oldY-1))
            for i in frange(1,math.ceil(yScale),step):
                tests.append((0,oldY-1-i))
        
        else:
            for i in frange(math.floor((p[1]-1)*yScale), math.ceil((p[1]+1)*yScale),step):
                tests.append((0,i))
    
    elif p[0] == newX - 1:
        if p[1] == 0:
            for i in frange(0,math.ceil(xScale),step):
                tests.append((oldX-1-i,0))
            for i in frange(1,math.ceil(yScale),step):
                tests.append((oldX-1,i))
        
        elif p[1] == newY - 1:
            for i in frange(0,math.ceil(xScale),step):
                tests.append((oldX-1-i,oldY-1))
            for i in frange(1,math.ceil(yScale),step):
                tests.append((oldX-1,oldY-1-i))
        
        else:
            for i in frange(math.floor((p[1]-1)*yScale), math.ceil((p[1]+1)*yScale),step):
                tests.append((oldX-1,i))

    else:
        if p[1] == 0:
            for i in frange(math.floor((p[0]-1)*xScale), math.ceil((p[0]+1)*xScale),step):
                tests.append((i,0))
        
        else:
            for i in frange(math.floor((p[0]-1)*xScale), math.ceil((p[0]+1)*xScale),step):
                tests.append((i,oldY-1))

    return tests


def findSatellite(imdata, fit):
    maxCount = -1
    p1 = (0,0)
    p2 = (0,0)

    minLen = max(imdata.shape[0],imdata.shape[1])/4
    
    for i in range(0, imdata.shape[1]-1):
        for j in range(1, imdata.shape[1]):
            count = searchLine(fit, (i,0), (j,imdata.shape[0]-1), minLen, median = True)
            if count>maxCount:
                maxCount = count
                p1 = (i,0)
                p2 = (j,imdata.shape[0]-1)

    
    for i in range(0, imdata.shape[0]-1):
        for j in range(1, imdata.shape[0]):
            count = searchLine(fit, (0,j), (imdata.shape[1]-1,i), minLen, median = True)
            if count>maxCount:
                maxCount = count
                p1 = (0,j)
                p2 = (imdata.shape[1]-1,i)

    for i in range(0, imdata.shape[0]-1):
        for j in range(0, imdata.shape[1]-1):
            count = searchLine(fit, (j,0), (0,i+1), minLen, median = True)
            if count>maxCount:
                maxCount = count
                p1 = (j,0)
                p2 = (0,i+1)

            count = searchLine(fit, (j+1,imdata.shape[0]-1), (0,i+1), minLen, median = True)
            if count>maxCount:
                maxCount = count
                p1 = (j+1,imdata.shape[0]-1)
                p2 = (0,i+1)

            count = searchLine(fit, (j,0), (imdata.shape[1]-1,i), minLen, median = True)
            if count>maxCount:
                maxCount = count
                p1 = (j,0)
                p2 = (imdata.shape[1]-1,i)

            count = searchLine(fit, (j+1,imdata.shape[0]-1), (imdata.shape[1]-1,i), minLen, median = True)
            if count>maxCount:
                maxCount = count
                p1 = (j+1,imdata.shape[0]-1)
                p2 = (imdata.shape[1]-1,i)

    return (p1,p2,fit)


def searchLine(fit, p1, p2, minLen, median = False):
    delX = p1[0]-p2[0]
    delY = p1[1]-p2[1]
    
    if abs(delX)>abs(delY):
        return searchLineX(fit, p1, p2, minLen, median)
    
    return searchLineY(fit, p1, p2, minLen, median)
    

def searchLineX(fit, p1, p2, minLen, median):
    slope = (p1[1]-p2[1])/(p1[0]-p2[0])
    if median:
        if p1[0] <= p2[0]:
            xs = np.arange(p1[0],p2[0]+1,1.0)
            ys = (xs-p1[0])*slope+p1[1]
            counts = [fit(xs[i],ys[i]) for i in range(0,xs.shape[0])]


        else:
            xs = np.arange(p2[0],p1[0]+1,1.0)
            ys = (xs-p2[0])*slope+p2[1]
            counts = [fit(xs[i],ys[i]) for i in range(0,xs.shape[0])]

        if len(counts)>minLen:
            return np.mean(counts) + np.median(counts)
        return 0

    else:
        if p1[0] <= p2[0]:
            xs = np.arange(p1[0],p2[0]+1,1.0)
            ys = (xs-p1[0])*slope+p1[1]
            counts = [fit(xs[i],ys[i]) for i in range(0,xs.shape[0])]


        else:
            xs = np.arange(p2[0],p1[0]+1,1.0)
            ys = (xs-p2[0])*slope+p2[1]
            counts = [fit(xs[i],ys[i]) for i in range(0,xs.shape[0])]

        if len(counts)>minLen:
            return np.mean(counts) 
        return 0

def searchLineY(fit, p1, p2, minLen, median):
    slope = (p1[0]-p2[0])/(p1[1]-p2[1])
    if median:
        if p1[1] <= p2[1]:
            ys = np.arange(p1[1],p2[1]+1,1.0)
            xs = (ys-p1[1])*slope+p1[0]
            counts = [fit(xs[i],ys[i]) for i in range(0,xs.shape[0])]


        else:
            ys = np.arange(p2[1],p1[1]+1,1.0)
            xs = (ys-p2[1])*slope+p2[0]
            counts = [fit(xs[i],ys[i]) for i in range(0,xs.shape[0])]

        if len(counts)>minLen:
            return np.mean(counts) + np.median(counts)
        return 0

    else:
        if p1[1] <= p2[1]:
            ys = np.arange(p1[1],p2[1]+1,1.0)
            xs = (ys-p1[1])*slope+p1[0]
            counts = [fit(xs[i],ys[i]) for i in range(0,xs.shape[0])]


        else:
            ys = np.arange(p2[1],p1[1]+1,1.0)
            xs = (ys-p2[1])*slope+p2[0]
            counts = [fit(xs[i],ys[i]) for i in range(0,xs.shape[0])]

        if len(counts)>minLen:
            return np.mean(counts) 
        return 0
