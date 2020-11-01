import numpy as np
import math

def gaussian(dists, brightness, sDev):
    perpDists = dists[1]
    nExp = -np.square(perpDists)/(2*(sDev**2))
    trail = brightness(dists[0])*np.power(math.e,nExp)

    return trail

def gaussian2(dists, b1, b2, b3, b4, b5, sDev):
    perpDists = dists[1]
    brightness = np.poly1d([b1, b2, b3, b4, b5])
    nExp = -np.square(perpDists)/(2*(sDev**2))
    trail = brightness(dists[0])*np.power(math.e,nExp)

    return trail

def getDistancesHorizontal(slope, intercept, xDim, yDim, offset = None):
    x = np.arange(0.0, xDim, 1.0)
    y = np.arange(0.0, yDim, 1.0)
    x,y = np.meshgrid(x,y)

    return findDistances(x,y,slope,intercept,offset)

def findDistances(x,y,slope,intercept,offset=None):
    delXProj = (x+slope*y-slope*intercept)/(1+slope**2)
    delYProj = slope*delXProj

    parDist = np.sqrt(np.square(delXProj)+np.square(delYProj))
    perpDist = (y-slope*x-intercept)/math.sqrt(1+slope**2)
    if offset!=None:
        return (parDist,perpDist-offset(parDist))
    return (parDist,perpDist)

def getDistancesVertical(x, xDim, yDim, offset=None):
    parDist = np.arange(0.0, yDim, 1.0)
    perpDist = np.arange(0.0, xDim, 1.0)-x
    perpDist, parDist = np.meshgrid(perpDist, parDist)
    
    if offset!=None:
        perpDist = perpDist-offset(parDist)
        
    return parDist,perpDist

def getXY(parDist,perpDist,slope,intercept,offset=None):
    scale = np.sqrt(slope**2+1)
    parXDist = 1/scale
    parYDist = slope/scale
    perpXDist = -slope/scale
    perpYDist = 1/scale
    
    if offset!=None:
        perpDist = perpDist+offset(parDist)
        
    x = parXDist*parDist+perpXDist*perpDist
    y = parYDist*parDist+perpYDist*perpDist+intercept
    
    return (x,y)
        
def generateTrailHorizontal(slope, intercept, brightness, sDev, xDim, yDim, offset = None):
    parDist, perpDist = getDistancesHorizontal(slope, intercept, xDim, yDim, offset)

    return gaussian((parDist, perpDist), brightness, sDev)

def generateTrailVertical(x, brightness, sDev, xDim, yDim, offset = None):
    parDist, perpDist = getDistancesVertical(x, xDim, yDim, offset)

    return gaussian((parDist, perpDist), brightness, sDev)
