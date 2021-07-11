import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import math
import os
import StarRemovalCode.PlateSolve as PS
import TrailRemovalCode.TrailFit as TF
import TrailRemovalCode.GenerateTrail as GT
from astroquery.vizier import Vizier
from astropy import coordinates
from astropy import units as u
from astropy import wcs
from astropy.table import Table
from astropy.nddata import NDData
from astropy.visualization import simple_norm
from astropy.coordinates import Angle
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
import photutils
from photutils.psf import DAOPhotPSFPhotometry as DAOP
from photutils.psf import IntegratedGaussianPRF as PRF
from photutils.psf import DAOGroup
from photutils.psf import BasicPSFPhotometry
from photutils.background import MMMBackground
from photutils.background import Background2D
from photutils.psf import extract_stars
from photutils import EPSFBuilder
from photutils import centroid_sources,centroid_com
import pickle


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

def showImageWithDetectedStars(img, title, starX, starY, starR, color, inverted):
    mx = img.max()
    mn = img.min()
    median = np.median(img)

    if inverted:
        plt.imshow(-img, cmap=cm.gray, vmin = -(mx+rescale*median)/(rescale+1), vmax = -(mn+rescale2*median)/(rescale2+1))
    else:
        plt.imshow(img, cmap=cm.gray, vmin = (mn+rescale2*median)/(rescale2+1), vmax = (mx+rescale*median)/(rescale+1))

    theta = np.arange(0,2*math.pi+math.pi/50,math.pi/50)
    for i in range(0,len(starR)):
        delX = starR[i]*np.cos(theta)
        delY = starR[i]*np.sin(theta)
        xPlot = starX[i]+delX
        yPlot = starY[i]+delY
            
        plt.plot(xPlot,yPlot,color+"-")

    plt.title(title)
    plt.show()

def showImageWithDetectedHotPixels(img, title, hotX, hotY, color, inverted):
    mx = img.max()
    mn = img.min()
    median = np.median(img)

    if inverted:
        plt.imshow(-img, cmap=cm.gray, vmin = -(mx+rescale*median)/(rescale+1), vmax = -(mn+rescale2*median)/(rescale2+1))
    else:
        plt.imshow(img, cmap=cm.gray, vmin = (mn+rescale2*median)/(rescale2+1), vmax = (mx+rescale*median)/(rescale+1))

    plt.plot(hotX,hotY,color+"o")

    plt.title(title)
    plt.show()

def findRadius(img, x, y, median, resolution = 50):
    r=1

    while True:
        counts = []
        for theta in range(0,resolution,1):
            boundX = int(round(x+r*math.cos(2*theta*math.pi/resolution)))
            boundY = int(round(y+r*math.sin(2*theta*math.pi/resolution)))

            try:
                counts.append(img[boundY][boundX])
            except:
                pass
        avg = np.median(counts)
        if avg <= median:
            return (r,avg)
        r+=1

def findRadiusGalaxy(img, x, y, median, resolution = 50):
    r=1

    thetas = np.array([theta for theta in range(0,resolution,1)])
    masks = []
    while True:
        counts = []
        for theta in thetas:
            boundX = int(round(x+r*math.cos(2*theta*math.pi/resolution)))
            boundY = int(round(y+r*math.sin(2*theta*math.pi/resolution)))

            try:
                counts.append(img[boundY][boundX])
            except:
                pass
        
        counts = np.array(counts)
        mask = counts>median
        if len(masks)!=0:
            for item in masks:
                mask = np.logical_or(mask,item)

        if len(masks)==2:
            masks = []
            thetas = thetas[mask]
        else:
            masks.append(mask)

        if thetas.shape[0] == 0:
            return r

        r+=1


def getStarCoords(img, pos, showProgress):
    if showProgress:
        print("Locating stars in the image...")

    v = Vizier(row_limit = 100000)#, columns = ['RA_ICRS', 'DE_ICRS','__Gmag_']

    coords = wcs.WCS(pos)
    center = coords.all_pix2world([(img.shape[0]/2,img.shape[1]/2)],0, ra_dec_order = True)[0]
    corner = coords.all_pix2world([(0,0)],0, ra_dec_order = True)[0]
    angle = np.sqrt((corner[0]-center[0])**2+(corner[1]-center[1])**2)

    c = coordinates.SkyCoord(center[0], center[1], unit=('deg', 'deg'), frame='icrs')
    starResult = v.query_region(c, radius=angle*u.deg, catalog = ["I/337/gaia","I/337/gaia2"])#["I/337/tgas","I/337/qso","I/345/"])
    galaxyNebulaResult = v.query_region(c, radius=angle*u.deg, catalog = ["VII/118/ngc2000"])

    print(starResult)

    starCoords = []
    mags = []
    for item in starResult:
        for i in range(0,item["RA_ICRS"].shape[0]):
            starCoords.append((item["RA_ICRS"][i],item["DE_ICRS"][i]))
            mags.append(item["__Gmag_"][i])

    starPix =  np.transpose(coords.all_world2pix(starCoords,0, ra_dec_order = True))
    starPixX = starPix[0]
    print(starPixX.shape[0])
    starPixY = starPix[1]

    showImageWithDetectedHotPixels(img, "Image with Stars", starPixX, starPixY, "b", True)

    mags = np.array(mags)

    galaxyCoords = []
    for item in galaxyNebulaResult:
        print(item)
        for i in range(0,item["RAB2000"].shape[0]):
            galaxyCoords.append((Angle(item["RAB2000"][i]+" hours").degree,Angle(item["DEB2000"][i]+" degrees").degree))

    galaxyPix =  np.transpose(coords.all_world2pix(galaxyCoords,0, ra_dec_order = True))
    if galaxyPix.shape[0] != 0:
        galaxyPixX = galaxyPix[0]
        galaxyPixY = galaxyPix[1]
    else:
        galaxyPixX = np.array([])
        galaxyPixY = np.array([])
    print(galaxyPixX)
    print(galaxyPixY)

    mask = (starPixX>0)

    mask = np.logical_and(mask, (starPixX<img.shape[1]))


    mask = np.logical_and(mask, (starPixY>0))
    mask = np.logical_and(mask, (starPixY<img.shape[0]))
    
    starPixX = starPixX[mask]
    starPixY = starPixY[mask]
    mags = mags[mask]

    order = np.argsort(mags)
    starPixX = starPixX[order]
    starPixY = starPixY[order]
    mags = mags[order]

    mask = (galaxyPixX>0)

    mask = np.logical_and(mask, (galaxyPixX<img.shape[1]))


    mask = np.logical_and(mask, (galaxyPixY>0))
    mask = np.logical_and(mask, (galaxyPixY<img.shape[0]))
    
    galaxyPixX = galaxyPixX[mask]
    galaxyPixY = galaxyPixY[mask]
    
    #magsS = np.sort(mags)
    #maxMag = magsS[3000]

    #mask = (mags<maxMag)
    #starPixX = starPixX[mask]
    #starPixY = starPixY[mask]
    #mags = mags[mask]
    
    if showProgress:
        print("Stars Located!\n")

    outFile = open("stars.st","wb")
    pickle.dump([starPixX,starPixY,mags,galaxyPixX,galaxyPixY],outFile)
    outFile.close()

    return (starPixX,starPixY,mags,galaxyPixX,galaxyPixY)

def removeGalaxies(img, gX, gY, showProgress, inverted = False, resolution = 100):
    background = np.median(img)#(np.median(img)+np.mean(img))/2
    
    showImage(img, "Image With Galaxies and Nebulae", inverted)

    if showProgress:
        print("Fitting Galaxy and Nebula Radii...")

    imgr = []
    for i in range(0,len(gX)):
        x= gX[i]
        y= gY[i]

        r = findRadiusGalaxy(img, x, y, background, resolution)
        imgr.append(r)

    if showProgress:
        print("Galaxy Radii Fitted!\n")
        showImageWithDetectedStars(img, "Image With Detected Galaxies and Nebulae", gX, gY, imgr, "r", inverted)

    if showProgress:
        print("Removing galaxies and nebulae from the image...")

    for i in range(0,len(gX)):
        x= gX[i]
        y= gY[i]
        xInt = int(x)
        yInt = int(y)

        r = imgr[i]

        for i in range(-r,r):
            h = round(math.sqrt(r**2-i**2))
            for j in range(-h,h):
                try:
                    img[yInt+i][xInt+j]=background
                except:
                    print("here")
                    pass

    if showProgress:
        print("Galaxies and Nebulae removed!\n")

        showImage(img, "Image With No Galaxies or Nebulae", inverted)

    return (img,imgr)

def removeStars(img, starX, starY, mags, gX, gY, gR, showProgress, inverted = False, resolution = 50):
    mask = np.ones(starX.shape[0],dtype=bool)
    for i in range(gX.shape[0]):
        mask = np.logical_and(mask, np.sqrt((starX-gX[i])**2+(starY-gY[i])**2)>gR[i])

    starPixX = starX[mask]
    starPixY = starY[mask]
    mags = mags[mask]

    background = np.mean(img)#(np.median(img)+np.mean(img))/2
    
    showImage(img, "Image With Stars", inverted)

    if showProgress:
        print("Fitting Star Radii...")

    imgr = []
    for i in range(0,len(starPixX)):
        x= starPixX[i]
        y= starPixY[i]

        r,avg = findRadius(img, x, y, background)
        imgr.append(r)

    if showProgress:
        print("Star Radii Fitted!\n")
        showImageWithDetectedStars(img, "Image With Detected Stars", starPixX, starPixY, imgr, "r", inverted)

    if showProgress:
        print("Removing stars from the image...")

    for i in range(0,len(starPixX)):
        x= starPixX[i]
        y= starPixY[i]
        xInt = int(x)
        yInt = int(y)

        r = imgr[i]

        for i in range(-r,r):
            h = round(math.sqrt(r**2-i**2))
            for j in range(-h,h):
                try:
                    img[yInt+i][xInt+j]=avg
                except:
                      pass

    if showProgress:
        print("Stars removed!\n")

        showImage(img, "Image With No Stars", inverted)

    return (img,starPixX,starPixY,mags,np.array(imgr))

def removeStarsPSF(img, starPixX, starPixY, starR, mags, gX, gY, gR, p1, p2, radius, showProgress, inverted = False):
    remove = []
    for i in range(0,len(starPixX)):
        if i%100 == 0:
            print(i)
        for j in range(i+1,len(starPixX)):
            if np.abs(starPixX[i]-starPixX[j])<2.5*radius and np.abs(starPixY[i]-starPixY[j])<2.5*radius:
                remove.append(i)
                remove.append(j)
    
    mask = np.ones(len(starPixX), dtype = bool)
    mask[np.array(remove)] = 0

    starPixX = starPixX[mask]
    starPixY = starPixY[mask]
    starR = starR[mask]
    mags = mags[mask]

    print(radius)
    mask = np.ones(starPixX.shape[0],dtype=bool)
    for i in range(gX.shape[0]):
        mask = np.logical_and(mask, np.sqrt((starPixX-gX[i])**2+(starPixY-gY[i])**2)>gR[i]+2*radius)

    starPixX = starPixX[mask]
    starPixY = starPixY[mask]
    starR = starR[mask]
    mags = mags[mask]

    mask = (2*radius<starPixX) & (starPixX<img.shape[1]-2*radius) & (2*radius<starPixY) & (starPixY<img.shape[0]-2*radius)
    starPixX = starPixX[mask]
    starPixY = starPixY[mask]
    starR = starR[mask]
    mags = mags[mask]

    if p1[0]==p2[0]:
        psfMask = np.abs(starPixX-p1[0])>2.5*radius
        trailMask = np.abs(starPixX-p1[0])<10
        psfX = starPixX[psfMask]
        psfY = starPixY[psfMask]
        psfR = starR[psfMask]
        psfMags = mags[psfMask]
        trailX = starPixX[trailMask]
        trailY = starPixY[trailMask]
        trailMags = mags[trailMask]
    else:
        m,b = TF.getLineParams(p1,p2)
        par,perp = GT.findDistances(starPixX,starPixY,m,b)
        psfMask = np.abs(perp)>2.5*radius
        trailMask = np.abs(perp)<10
        psfX = starPixX[psfMask]
        psfY = starPixY[psfMask]
        psfR = starR[psfMask]
        psfMags = mags[psfMask]
        trailX = starPixX[trailMask]
        trailY = starPixY[trailMask]
        trailMags = mags[trailMask]

    mask = psfR>10

    psfX = psfX[mask]
    psfY = psfY[mask]
    psfMags = psfMags[mask]
    print(psfX.shape[0])

    showImageWithDetectedHotPixels(img, "Image with Star Picks", psfX, psfY, "b", inverted)

    psfX,psfY = centroid_sources(img, psfX, psfY, box_size=21, centroid_func=centroid_com)

    showImageWithDetectedHotPixels(img, "Image with Star Picks Centroid", psfX, psfY, "b", inverted)

    stars = Table()
    stars["x"] = psfX
    stars["y"] = psfY

    #stars2 = Table()
    #stars2["x_0"] = psfX
    #stars2["y_0"] = psfY

    stars = extract_stars(NDData(data=img), stars, size = 2*int(2*radius)+1)

    nrows = 5
    ncols = 5
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
    ax = ax.ravel()
    for i in range(nrows*ncols):
        norm = simple_norm(stars[i], 'log', percent=99.)
        ax[i].imshow(stars[i], norm=norm, origin='lower', cmap='viridis')

    plt.show()

    #sigma_psf = 2.0
    #daogroup = DAOGroup(2.0)
    #mmm_bkg = MMMBackground()
    #fitter = LevMarLSQFitter()
    #psf_model = PRF(sigma=sigma_psf)
    #photometry = BasicPSFPhotometry(group_maker=daogroup,bkg_estimator=mmm_bkg,psf_model=psf_model,fitter=LevMarLSQFitter(),fitshape=(11,11))
    #result_tab = photometry(image=img, init_guesses=stars2)
    #residual_image = photometry.get_residual_image()

    epsfBuilder = EPSFBuilder(oversampling=2, maxiters=20, smoothing_kernel = "quadratic")
    epsf, starFits = epsfBuilder(stars)

    #showImage(residual_image,"Image No Stars",True)

    norm = simple_norm(epsf.data, 'log', percent=99.)
    plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
    plt.colorbar()
    plt.show()

    nrows = 5
    ncols = 5
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
    ax = ax.ravel()
    for i in range(nrows*ncols):
        norm = simple_norm(starFits[i], 'log', percent=99.)
        ax[i].imshow(starFits[i], norm=norm, origin='lower', cmap='viridis')

    plt.show()

    nrows = 5
    ncols = 5
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20), squeeze=True)
    ax = ax.ravel()
    for i in range(nrows*ncols):
        norm = simple_norm(stars[i].data-starFits[i].data, 'log', percent=99.)
        ax[i].imshow(stars[i].data-starFits[i].data, norm=norm, origin='lower', cmap='viridis')

    plt.show()

    showImage(finalImg, "Image With No Stars", inverted)


def removeOutlierPixels(img, showProgress, inverted = False):
    if showProgress:
        print("Removing Outlier Pixels...")
    
    sDev = np.std(img)
    mean = np.mean(img)
    zScore = np.abs(img-mean)/sDev
    mask = (zScore > 4)

    outlierCandidatesY,outlierCandidatesX = mask.nonzero()

    if showProgress:
        showImageWithDetectedHotPixels(img, "Image with Outlier Candidates", outlierCandidatesX, outlierCandidatesY, "b", inverted)

    confirmedOutliersX = []
    confirmedOutliersY = []
    boundCheck = 2
    for i in range(0,len(outlierCandidatesX)):
        xs = []
        ys = []

        xMin = max(0,outlierCandidatesY[i]-boundCheck)
        xMax = min(img.shape[0],outlierCandidatesY[i]+boundCheck+1)
        yMin = max(0,outlierCandidatesX[i]-boundCheck)
        yMax = min(img.shape[1],outlierCandidatesX[i]+boundCheck+1)

        box = img[xMin:xMax,yMin:yMax]


        bMean = np.mean(box)
        bSDev = np.std(box)

        bZScore = np.abs(img[outlierCandidatesY[i],outlierCandidatesX[i]]-bMean)/bSDev

        if (bZScore >= 3):
            img[outlierCandidatesY[i],outlierCandidatesX[i]] = bMean
            if showProgress:
                confirmedOutliersX.append(outlierCandidatesX[i])
                confirmedOutliersY.append(outlierCandidatesY[i])



    if showProgress:
        showImageWithDetectedHotPixels(img, "Image with Identified Outliers", confirmedOutliersX, confirmedOutliersY, "b", inverted)
        print("Outlier Pixels Removed!\n\n")
        showImage(img, "Image With No Stars Or Outliers", inverted)

    return (img, confirmedOutliersX, confirmedOutliersY)


def findAndRemoveStars(img, imgPath, showProgress = True, inverted = False):
    showImage(img, "Image With Stars", inverted)
    #pos = PS.getStars(imgPath, showProgress)
    #starPixX,starPixY,mags,gX,gY = getStarCoords(img, pos, showProgress)
    fileIn = open("stars.st","rb")
    starPixX,starPixY,mags,gX,gY = pickle.load(fileIn)
    fileIn.close()
    #print(gX)
    #print(gY)
    imgNoStars, gR = removeGalaxies(img, gX, gY, showProgress, inverted = inverted, resolution = 50)
    imgNoStars, starPixX, starPixY, mags, starR = removeStars(imgNoStars, starPixX, starPixY, mags, gX, gY, gR, showProgress, inverted = inverted)
    imgNoStarsNoOutliers, hotX, hotY = removeOutlierPixels(imgNoStars, showProgress, inverted = inverted)

    return (imgNoStarsNoOutliers, gX, gY, gR, starPixX, starPixY, starR, mags, hotX, hotY)
