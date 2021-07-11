import numpy as np
import math
from numpy import save
import TrailRemovalCode.GenerateTrail as GT
import TrailRemovalCode.RemoveTrail as RT
import TrailRemovalCode.FourierFit as FF
from scipy.optimize import curve_fit as cf

from photutils.background import MMMBackground

from matplotlib import cm
from matplotlib import pyplot as plt

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

def showImageWithNearbyStars(img, title, starX, starY, starR, inverted):
    mx = img.max()
    mn = img.min()
    median = np.median(img)

    if inverted:
        plt.imshow(-img, cmap=cm.gray, vmin = -(mx+rescale*median)/(rescale+1), vmax = -(mn+rescale2*median)/(rescale2+1))
    else:
        plt.imshow(img, cmap=cm.gray, vmin = (mn+rescale2*median)/(rescale2+1), vmax = (mx+rescale*median)/(rescale+1))
    
    plt.title(title)

    for i in range(0,len(starX)):
        theta = np.arange(0,2*math.pi+math.pi/50,math.pi/50)
        delX = starR[i]*np.cos(theta)
        delY = starR[i]*np.sin(theta)
        xPlot = starX[i]+delX
        yPlot = starY[i]+delY
                
        plt.plot(xPlot,yPlot,"b-")

    plt.show()

def getLineParams(p1, p2):
    m = (p2[1]-p1[1])/(p2[0]-p1[0])
    b = p1[1]-m*p1[0]
    return (m,b)

def getParLine(img, p1, p2, m, b):
    parDists = GT.findDistances(np.array([p1[0],p2[0]]),np.array([p1[1],p2[1]]),m,b)[0]

    minD = min(parDists)
    maxD = max(parDists)
    step = (maxD-minD)/max(img.shape)

    return np.arange(minD,maxD,step)

def getParLineXYBand(perpLine, img, p1, p2, m, b, offset=None):
    parLine = getParLine(img, p1, p2, m, b)

    parLine,perpLine = np.meshgrid(parLine,perpLine)
    return GT.getXY(parLine,perpLine,m,b,offset)

def getStarsHotPixelsAlongTrail(m, b, starX, starY, starR, gX, gY, gR, hotX, hotY, dist):
    starX = np.array(starX)
    starY = np.array(starY)
    starR = np.array(starR)
    gX = np.array(gX)
    gY = np.array(gY)
    gR = np.array(gR)
    hotX = np.array(hotX)
    hotY = np.array(hotY)

    starPar,starPerp = GT.findDistances(starX,starY,m,b)
    gPar,gPerp = GT.findDistances(gX,gY,m,b)
    hotPar,hotPerp = GT.findDistances(hotX,hotY,m,b)

    starMask = (np.abs(np.abs(starPerp)-np.abs(starR))<=dist)
    gMask = (np.abs(np.abs(gPerp)-np.abs(gR))<=dist)
    hotMask = (np.abs(hotPerp)<=dist)

    return (starX[starMask], starY[starMask], starR[starMask], gX[gMask], gY[gMask], gR[gMask], hotX[hotMask], hotY[hotMask])

def getTrailBand(parDist, perpDist, img, dist):
    parDistF = parDist.flatten()
    perpDistF = perpDist.flatten()
    imgF = img.flatten()
    mask = (np.abs(perpDistF) <= dist)
    perpDistF = perpDistF[mask]
    parDistF = parDistF[mask]
    imgF = imgF[mask]

    return (parDistF,perpDistF,imgF)


def getNonTrailBand(parDist, perpDist, img, dist1, dist2):
    parDistF = parDist.flatten()
    perpDistF = perpDist.flatten()
    imgF = img.flatten()
    mask = (np.abs(perpDistF) >= dist1)
    mask = np.logical_and(mask,(np.abs(perpDistF) <= dist2))

    perpDistF = perpDistF[mask]
    parDistF = parDistF[mask]
    imgF = imgF[mask]

    return (parDistF,perpDistF,imgF)


def interpolateSignal(x, y, mask):
    print("Beggining Trail Interpolation")
    prevI = -1

    for i in range(0,len(x)):
        if mask[i]:
            if i != prevI+1:
                if prevI == -1:
                    y[0:i] = y[i]
                else:
                    #print("Previous I: "+str(y[prevI]))
                    #print("Current I: "+str(y[i]))
                    m = (y[i]-y[prevI])/(x[i]-x[prevI])
                    b = y[i]-m*x[i]
                    y[prevI:i] = m*x[prevI:i]+b

            prevI = i

    print(y[prevI])
    print(y[prevI:])

    y[prevI:] = y[prevI]


def getValueMask(img, x, y, starX, starY, starR, gX, gY, gR, hotX, hotY, minDist):
    mask = np.ones(len(x),dtype = bool)
    for i in range(0,len(starX)):
        mask = np.logical_and(mask, (np.sqrt((x-starX[i])**2+(y-starY[i])**2) > starR[i]+minDist))

    for i in range(0,len(gX)):
        mask = np.logical_and(mask, (np.sqrt((x-gX[i])**2+(y-gY[i])**2) > gR[i]+minDist))

    for i in range(0,len(hotX)):
        mask = np.logical_and(mask, (np.sqrt((x-hotX[i])**2+(y-hotY[i])**2) > 1+minDist))

    mask = np.logical_and(mask, x>=0)

    mask = np.logical_and(mask, y>=0)

    mask = np.logical_and(mask, x<=img.shape[1])

    mask = np.logical_and(mask, y<=img.shape[0])

    return mask



def findOffsetX(p1, p2, m, b, starX, starY, starR, hotX, hotY, img, PC):
    perpLine = np.arange(-4,4.1,0.1)
    parLine = getParLine(img, p1, p2, m, b)
    x,y = getParLineXYBand(perpLine, img, p1, p2, m, b)

    pars = []
    offsets = []
    for i in range(0,x.shape[1]):
        bright = [PC(x[j,i],y[j,i]) for j in range(0,x.shape[0])]
        off = [perpLine[j] for j in range(0,x.shape[0])]
        maxOff = off[np.argmax(bright)]

        offsets.append(maxOff)
        pars.append(parLine[i])

    offsets = np.array(offsets)

    parLine2 = getParLine(img, p1, p2, m, b)
    trailX,trailY = GT.getXY(parLine2, offsets, m, b)

    fig,axes = plt.subplots(2,1,True,True)
    axes[0].plot(pars,offsets)
    axes[0].set_title("Wobble with Stars")

    mask = getValueMask(img, trailX, trailY, starX, starY, starR, hotX, hotY, 1)
    interpolateSignal(parLine2, offsets, mask)

    axes[1].plot(pars,offsets)
    axes[1].set_title("Wobble without Stars")
    plt.show()


    pars = np.array(pars)
    offset = FF.FourierFit(pars, offsets, int(pars.shape[0]*0.87))#0.92))
    
    if True:
        plt.plot(pars,offsets,label="Wobble")
        plt.plot(pars,offset(pars), label="Wobble Fit")
        plt.legend()
        plt.show()

    x2,y2 = GT.getXY(pars, offset(pars), m, b)

    if True:
        mx = img.max()
        mn = img.min()
        median = np.median(img)

        plt.imshow(-img, cmap=cm.gray, vmin = -(mx+rescale*median)/(rescale+1), vmax = -(mn+rescale2*median)/(rescale2+1))
        plt.title("Trail With Wobble Fit")

        plt.plot(x2,y2)
        
        plt.show()

    return offset


def findOffsetX2(p1, p2, m, b, starX, starY, starR, gX, gY, gR, hotX, hotY, img, PC, perpRangeSize, step, gaussRange, sDev):
    perpLine = np.arange(-perpRangeSize-gaussRange, perpRangeSize+gaussRange+step, step)
    parLine = getParLine(img, p1, p2, m, b)
    x,y = getParLineXYBand(perpLine, img, p1, p2, m, b)

    perpRange = np.arange(-gaussRange,gaussRange+step,step)
    gaussian = gauss(perpRange,sDev)

    plt.plot(perpRange,gaussian)
    plt.show()

    startI = int(np.round(gaussRange/step))
    endI = int(np.round(len(perpLine)-1-gaussRange/step))

    pars = []
    offsets = []
    for i in range(0, x.shape[1]):
        bright = [PC(x[j,i],y[j,i]) for j in range(0,x.shape[0])]
        gaussBright = []

        for j in range(startI,endI):
            brightList = []
            for k in range(-startI,startI):
                brightList.append(bright[j+k]*gaussian[startI+k])
            gaussBright.append(np.mean(brightList))

        off = [perpLine[j] for j in range(startI,endI)]
        maxOff = off[np.argmax(gaussBright)]

        offsets.append(maxOff)
        pars.append(parLine[i])

    offsets = np.array(offsets)

    parLine2 = getParLine(img, p1, p2, m, b)
    trailX,trailY = GT.getXY(parLine2, offsets, m, b)

    fig,axes = plt.subplots(2,1,True,True)
    axes[0].plot(pars,offsets)
    axes[0].set_title("Wobble with Stars")

    mask = getValueMask(img, trailX, trailY, starX, starY, starR, gX, gY, gR, hotX, hotY, 1)
    interpolateSignal(parLine2, offsets, mask)

    axes[1].plot(pars,offsets)
    axes[1].set_title("Wobble without Stars")
    plt.show()

    pars = np.array(pars)
    offset = FF.FourierFit(pars, offsets, int(pars.shape[0]*0.95))
    
    if True:
        plt.plot(pars,offsets,label="Wobble")
        plt.plot(pars,offset(pars), label="Wobble Fit")
        plt.legend()
        plt.title("Wobble Fit")
        plt.show()

    x2,y2 = GT.getXY(pars, offset(pars), m, b)

    if True:
        mx = img.max()
        mn = img.min()
        median = np.median(img)

        plt.imshow(-img, cmap=cm.gray, vmin = -(mx+rescale*median)/(rescale+1), vmax = -(mn+rescale2*median)/(rescale2+1))
        plt.title("Trail With Wobble Fit")

        plt.plot(x2,y2)
        plt.plot(trailX,trailY)
        plt.plot([p1[0],p2[0]],[p1[1],p2[1]])
        
        plt.show()

    return offset
    

def fitBrightnessPolyX(img, p1, p2, m, b, offset, PC, starX, starY, starR, gX, gY, gR, hotX, hotY):
    parDists = getParLine(img, p1, p2, m, b)
    perpDists = parDists*0
    x,y = GT.getXY(parDists,perpDists,m,b,offset)

    #weights = ((dists-np.mean(dists))/max(dists))**2+1
        
    brightness = np.array([PC(x[i],y[i]) for i in range(0,len(x))])

    fig,axes = plt.subplots(2,1,True,True)
    axes[0].plot(parDists,brightness)
    axes[0].set_title("Britghtness with Stars")

    mask = getValueMask(img, x, y, starX, starY, starR, gX, gY, gR, hotX, hotY, 1)
    interpolateSignal(parDists, brightness, mask)

    axes[1].plot(parDists,brightness)
    axes[1].set_title("Brightness without Stars")
    plt.show()

    bPol = np.poly1d(np.polyfit(parDists,brightness,20))#,weights = weights))
    bFour = FF.FourierFit(parDists, brightness, int(parDists.shape[0]*0.98))#0.92))

    plt.plot(parDists, brightness, label = "Brightness")
    plt.plot(parDists, bPol(parDists), label = "Brightness Fit")
    #plt.plot(parDists, bFour(parDists), label = "Brightness Fit")
    plt.title("Brightness Fit")
    plt.legend()
    plt.show()

    return bFour

def fitBrightnessPolyX2(img, p1, p2, m, b, offset, PC, starX, starY, starR, hotX, hotY, gaussRange, step, sDev):
    parLine = getParLine(img, p1, p2, m, b)
    perpLine = np.arange(-gaussRange,gaussRange+step,step)
    x,y = getParLineXYBand(perpLine, img, p1, p2, m, b, offset)

    gaussian = gauss(perpLine,sDev)

    brightness = []
    for i in range(0,x.shape[1]):
        brightness.append(np.mean([PC(x[j,i],y[j,i])*gaussian[j] for j in range(0,len(x))]))

    brightness = np.array(brightness)

    trailX,trailY = GT.getXY(parLine, parLine*0, m, b, offset)
    mask = getValueMask(img, trailX, trailY, starX, starY, starR, hotX, hotY, gaussRange)

    interpolateSignal(parLine, brightness, mask)


    bPol = np.poly1d(np.polyfit(parLine,brightness,5))#,weights = weights))
    #bFour = FF.FourierFit(parLine, brightness, int(parLine.shape[0]*0.98))#0.92))

    plt.plot(parLine, brightness)
    plt.plot(parLine, bPol(parLine))
    #plt.plot(parLine, bFour(parLine))
    plt.title("Brightness Fit")
    plt.show()

    middleIndex = int(np.round(gaussRange/step))
    brightness2 = np.array([(PC(x[middleIndex,i],y[middleIndex,i])/bPol(parLine[i])) for i in range(0,len(x))])

    bPol2 = np.poly1d(np.polyfit(parLine,brightness2,0))#,weights = weights))
    #bFour = FF.FourierFit(parLine, brightness, int(parLine.shape[0]*0.98))#0.92))


    plt.plot(parLine, brightness2)
    plt.plot(parLine, bPol2(parLine))
    #plt.plot(parLine, bFour(parLine))
    plt.title("Brightness Fit")
    plt.show()

    return bPol

def estimateBackground(img, m, b, offset, starX, starY, starR, gX, gY, gR, hotX, hotY):
    parDist,perpDist = GT.getDistancesHorizontal(m, b, img.shape[1], img.shape[0], offset)
    par,perp,band = getNonTrailBand(parDist, perpDist, img, 10, 20)
    x,y = GT.getXY(par,perp,m,b,offset)

    mask = getValueMask(img, x, y, starX, starY, starR, gX, gY, gR, hotX, hotY, 1)
    par = par[mask]
    band = band[mask]

    mean = np.mean(band)

    plt.plot(par,band,"bo",label = "Background Pixel Brightnesses")
    plt.plot([par[0],par[-1]],[mean, mean],label = "Average Background Brightness: "+str(np.round(mean,decimals = 3)))
    plt.title("Background Estimate")
    plt.legend()
    plt.show()

    print(mean)
    return mean


def fitGaussianX(img, m, b, bright, offset, starX, starY, starR, gX, gY, gR, hotX, hotY):
    parDist,perpDist = GT.getDistancesHorizontal(m, b, img.shape[1], img.shape[0], offset)
    parF,perpF,imgF = getTrailBand(parDist, perpDist, img, 6)
    x,y = GT.getXY(parF,perpF,m,b,offset)
    
    trailScaled = imgF/bright(parF)

    fig,axes = plt.subplots(1,2,sharex = True, sharey = True)

    axes[0].plot(perpF,trailScaled,"bo")
    axes[0].set_title("Trail Gaussian With Stars")

    mask = getValueMask(img, x, y, starX, starY, starR, gX, gY, gR, hotX, hotY, 1)
    perpF = perpF[mask]
    trailScaled = trailScaled[mask]
    

    axes[1].plot(perpF,trailScaled,"bo")
    axes[1].set_title("Trail Gaussian Without Stars")

    plt.show()


    sDev = cf(gauss, perpF, trailScaled)[0][0]


    x = np.arange(-5,5,0.1)
    
    plt.plot(perpF,trailScaled,"bo",label = "Pixel Brightnesses")
    plt.plot(x,gauss(x,sDev), label = "Gaussian Fit")
    plt.title("Gaussian Fit")
    plt.show()

    return sDev

def fitBrightnessSDev(img, m, b, offset, starX, starY, starR, gX, gY, gR, hotX, hotY, trailSize):
    parDist,perpDist = GT.getDistancesHorizontal(m, b, img.shape[1], img.shape[0], offset)
    parDistF,perpDistF,imgF = getTrailBand(parDist, perpDist, img, trailSize)
    x,y = GT.getXY(parDistF,perpDistF,m,b,offset)

    mask = getValueMask(img, x, y, starX, starY, starR, gX, gY, gR, hotX, hotY, 1)
    par = parDistF[mask]
    perp = perpDistF[mask]
    band = imgF[mask]

    b1, b2, b3, b4, b5, sDev = cf(GT.gaussian2, [par, perp], band, bounds = ([-math.inf,-math.inf,-math.inf,-math.inf,-math.inf,0],[math.inf,math.inf,math.inf,math.inf,math.inf,trailSize]))[0]

    return [np.poly1d([b1, b2, b3, b4, b5]), sDev]
    

def findFitHorizontal(img, p1, p2, starX, starY, starR, gX, gY, gR, hotX, hotY, trailSize, PC):
    m, b = getLineParams(p1, p2)

    starX, starY, starR, gX, gY, gR, hotX, hotY = getStarsHotPixelsAlongTrail(m, b, starX, starY, starR, gX, gY, gR, hotX, hotY, 25)
    showImageWithNearbyStars(img, "Stars Detected Nearby", starX, starY, starR, True)
    
    offset = findOffsetX2(p1, p2, m, b, starX, starY, starR, gX, gY, gR, hotX, hotY, img, PC, 5, 0.1, 2, 2)
    
    back = estimateBackground(img, m, b, offset, starX, starY, starR, gX, gY, gR, hotX, hotY)
    img = img-back

    #brightness, sDev = fitBrightnessSDev(img, m, b, offset, starX, starY, starR, gX, gY, gR, hotX, hotY, trailSize)

    brightness = fitBrightnessPolyX(img, p1, p2, m, b, offset, PC, starX, starY, starR, gX, gY, gR, hotX, hotY)#, 2, 0.1, 2)
    
    sDev = fitGaussianX(img, m, b, brightness, offset, starX, starY, starR, gX, gY, gR, hotX, hotY)

    return (sDev,brightness,offset)
    

def findFitVertical(img, x, trailSize):
    parDist, perpDist = GT.getDistancesVertical(x, img.shape[1], img.shape[0])

    median = np.median(img)

    bkg = MMMBackground()
    back = bkg(img)

    return fit(img-median, perpDist, parDist, trailSize)

def gauss(x, sDev):
    return np.power(math.e,-np.square(x)/(2*(sDev**2)))
