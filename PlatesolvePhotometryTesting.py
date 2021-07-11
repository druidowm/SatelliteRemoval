import os
import StarRemovalCode.PlateSolve as PS
import TrailRemovalCode.GenerateTrail as GT
import numpy as np
import astropy as ap
from astropy import coordinates
from astropy import units as u
from astropy import wcs
from astropy.table import Table
from astropy.nddata import NDData
from astropy.visualization import simple_norm
from astropy.coordinates import Angle
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.io import fits
from astroquery.vizier import Vizier
from matplotlib import cm
from matplotlib import pyplot as plt
from photutils import Background2D
from photutils import CircularAperture as CAp
from photutils import CircularAnnulus as CAn
from photutils.aperture import aperture_photometry as ap
from photutils.utils import calc_total_error as cte
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

def getLineParams(p1, p2):
    m = (p2[1]-p1[1])/(p2[0]-p1[0])
    b = p1[1]-m*p1[0]
    return (m,b)

def loadImages(trailImg, noTrailImg, baselineImgs):
    trail = fits.open(trailImg)[0].data.astype(np.float32)

    noTrail = fits.open(noTrailImg)[0].data.astype(np.float32)
    
    images = os.scandir(baselineImgs)

    baselineImages = []
    baselineNames = []

    for image in images:
        if image.name != ".DS_Store" and not image.is_dir():
            baselineImages.append(fits.open(baselineImgs+"/"+image.name)[0].data.astype(np.float32))
            baselineNames.append(baselineImgs+"/"+image.name)

    showImage(trail,"",True)
    showImage(noTrail,"",True)
    for image in baselineImages:
        showImage(image,"",True)
    return (trail, noTrail, baselineImages, baselineNames)

def getPos(trailImg, noTrailImg, baseslineImages, baselineImgs):
    trailPos = PS.getStars(trailImg, False)
    noTrailPos = PS.getStars(noTrailImg, False)

    baselinePos = []

    for img in baselineImgs:
        baselinePos.append(PS.getStars(img, False))

    trailCoords = wcs.WCS(trailPos)
    noTrailCoords = wcs.WCS(noTrailPos)

    baselineCoords = []
    for pos in baselinePos:
        baselineCoords.append(wcs.WCS(pos))
    
    return (trailCoords,noTrailCoords,baselineCoords)

def getStars(trail, noTrail, baselineImages, trailCoords, noTrailCoords, baselineCoords, p1, p2):
    v = Vizier(row_limit = 100000)

    center = trailCoords.all_pix2world([(trail.shape[0]/2,trail.shape[1]/2)],0, ra_dec_order = True)[0]
    corner = trailCoords.all_pix2world([(0,0)],0, ra_dec_order = True)[0]
    angle = np.sqrt((corner[0]-center[0])**2+(corner[1]-center[1])**2)

    c = coordinates.SkyCoord(center[0], center[1], unit=('deg', 'deg'), frame='icrs')
    starResult = v.query_region(c, radius=angle*u.deg, catalog = ["I/337/gaia"])

    starCoords = []
    mags = []
    for item in starResult:
        for i in range(0,item["RA_ICRS"].shape[0]):
            starCoords.append((item["RA_ICRS"][i],item["DE_ICRS"][i]))
            mags.append(item["__Gmag_"][i])

    starPix =  np.transpose(trailCoords.all_world2pix(starCoords,0, ra_dec_order = True))
    starPixX = starPix[0]
    starPixY = starPix[1]
    mags = np.array(mags)

    mask = (starPixX>50)

    mask = np.logical_and(mask, (starPixX<trail.shape[1]-50))


    mask = np.logical_and(mask, (starPixY>50))
    mask = np.logical_and(mask, (starPixY<trail.shape[0]-50))
    mask = np.logical_and(mask, (starPixX+starPixY<6000))
    
    mags = mags[mask]

    newStarCoords = []
    for i in range(len(starCoords)):
        if mask[i]:
            newStarCoords.append(starCoords[i])
    starCoords = newStarCoords


    starPix =  np.transpose(noTrailCoords.all_world2pix(starCoords,0, ra_dec_order = True))
    starPixX = starPix[0]
    starPixY = starPix[1]
    mags = np.array(mags)

    mask = (starPixX>50)

    mask = np.logical_and(mask, (starPixX<trail.shape[1]-50))


    mask = np.logical_and(mask, (starPixY>50))
    mask = np.logical_and(mask, (starPixY<trail.shape[0]-50))
    
    mags = mags[mask]

    newStarCoords = []
    for i in range(len(starCoords)):
        if mask[i]:
            newStarCoords.append(starCoords[i])
    starCoords = newStarCoords

    for coords in baselineCoords:
        starPix =  np.transpose(coords.all_world2pix(starCoords,0, ra_dec_order = True))
        starPixX = starPix[0]
        starPixY = starPix[1]
        mags = np.array(mags)

        mask = (starPixX>50)

        mask = np.logical_and(mask, (starPixX<trail.shape[1]-50))


        mask = np.logical_and(mask, (starPixY>50))
        mask = np.logical_and(mask, (starPixY<trail.shape[0]-50))
        
        mags = mags[mask]

        newStarCoords = []
        for i in range(len(starCoords)):
            if mask[i]:
                newStarCoords.append(starCoords[i])
        starCoords = newStarCoords

    trailStarPix =  np.transpose(trailCoords.all_world2pix(starCoords,0, ra_dec_order = True))
    trailStarPixX = trailStarPix[0]
    trailStarPixY = trailStarPix[1]

    showImageWithDetectedHotPixels(trail, "", trailStarPixX, trailStarPixY, "b", True)


    noTrailStarPix =  np.transpose(noTrailCoords.all_world2pix(starCoords,0, ra_dec_order = True))
    noTrailStarPixX = noTrailStarPix[0]
    noTrailStarPixY = noTrailStarPix[1]

    showImageWithDetectedHotPixels(noTrail, "", noTrailStarPixX, noTrailStarPixY, "b", True)

    baselineStarPixX = []
    baselineStarPixY = []
    for i in range(len(baselineCoords)):
        baselinePix =  np.transpose(baselineCoords[i].all_world2pix(starCoords,0, ra_dec_order = True))
        baselineStarPixX.append(baselinePix[0])
        baselineStarPixY.append(baselinePix[1])
        showImageWithDetectedHotPixels(baselineImages[i], "", baselineStarPixX[-1], baselineStarPixY[-1], "b", True)

    order = np.argsort(mags)
    trailStarPixX = trailStarPixX[order]
    trailStarPixY = trailStarPixY[order]
    noTrailStarPixX = noTrailStarPixX[order]
    noTrailStarPixY = noTrailStarPixY[order]
    for i in range(len(baselineStarPixX)):
        baselineStarPixX[i] = baselineStarPixX[i][order]
        baselineStarPixY[i] = baselineStarPixY[i][order]
    mags = mags[order]

    if p1[0]==p2[0]:
        perp = trailStarPixX-p1[0]
    else:
        m,b = getLineParams(p1,p2)
        par,perp = GT.findDistances(trailStarPixX,trailStarPixY,m,b)

    mask = (np.abs(perp)<=20)
    antimask = (np.abs(perp)>20)

    baselineStarTrailPixX = []
    baselineStarTrailPixY = []

    trailStarTrailPixX = trailStarPixX[mask][:12]
    trailStarTrailPixY = trailStarPixY[mask][:12]

    showImageWithDetectedHotPixels(trail, "", trailStarTrailPixX, trailStarTrailPixY, "b", True)

    noTrailStarTrailPixX = noTrailStarPixX[mask][:12]
    noTrailStarTrailPixY = noTrailStarPixY[mask][:12]

    showImageWithDetectedHotPixels(noTrail, "", noTrailStarTrailPixX, noTrailStarTrailPixY, "b", True)

    for i in range(len(baselineStarPixX)):
        baselineStarTrailPixX.append(baselineStarPixX[i][mask][:12])
        baselineStarTrailPixY.append(baselineStarPixY[i][mask][:12])
        showImageWithDetectedHotPixels(baselineImages[i], "", baselineStarTrailPixX[-1], baselineStarTrailPixY[-1], "b", True)

    baselineStarBackPixX = []
    baselineStarBackPixY = []

    trailStarBackPixX = trailStarPixX[antimask][:200]
    trailStarBackPixY = trailStarPixY[antimask][:200]

    showImageWithDetectedHotPixels(trail, "", trailStarBackPixX, trailStarBackPixY, "b", True)

    for i in range(len(baselineStarPixX)):
        baselineStarBackPixX.append(baselineStarPixX[i][antimask][:200])
        baselineStarBackPixY.append(baselineStarPixY[i][antimask][:200])
        showImageWithDetectedHotPixels(baselineImages[i], "", baselineStarBackPixX[-1], baselineStarBackPixY[-1], "b", True)
    
    return (trailStarTrailPixX,trailStarTrailPixY,noTrailStarTrailPixX,noTrailStarTrailPixY,baselineStarTrailPixX,baselineStarTrailPixY,trailStarBackPixX,trailStarBackPixY,baselineStarBackPixX,baselineStarBackPixY)

def getAperturePhotometry(img, starX, starY):
    stars = [(starX[i],starY[i]) for i in range(starX.shape[0])]
    aperture = CAp(stars, r=8)

    aperture_masks = aperture.to_mask()
    star_apertures = [aperture_masks[i].multiply(img) for i in range(0,len(aperture_masks))]



    #annulus_masks = annulus.to_mask(method = 'center')
    #star_annulus = [annulus_masks[i].multiply(data2) for i in range(0,len(annulus_masks))]

    #for i in range(0,len(star_apertures)):
    #    plt.imshow(star_apertures[i])
    #    plt.show()

    phot_table = ap(img, aperture)

    return phot_table["aperture_sum"]

def main(trailImg, noTrailImg, baselineImgs, p1, p2, plateSolve = False):
    trail, noTrail, baselineImages, baselineNames = loadImages(trailImg, noTrailImg, baselineImgs)

    if plateSolve:
        trailCoords, noTrailCoords, baselineCoords = getPos(trailImg, noTrailImg, baselineImages, baselineNames)
        saveFile = open("starsPhotometry.st","wb")
        pickle.dump([trailCoords, noTrailCoords, baselineCoords], saveFile)
        saveFile.close()
    else:
        saveFile = open("starsPhotometry.st","rb")
        trailCoords, noTrailCoords, baselineCoords = pickle.load(saveFile)
        saveFile.close()

    TSTX,TSTY,NTSTX,NTSTY,BSTX,BSTY,TSBX,TSBY,BSBX,BSBY = getStars(trail, noTrail, baselineImages, trailCoords, noTrailCoords, baselineCoords, p1, p2)

    TSTPhot = getAperturePhotometry(trail, TSTX, TSTY)
    NTSTPhot = getAperturePhotometry(noTrail, NTSTX, NTSTY)
    BSTPhot = [getAperturePhotometry(baselineImages[i], BSTX[i], BSTY[i]) for i in range(len(baselineImages))]
    TSBPhot = getAperturePhotometry(trail, TSBX, TSBY)
    BSBPhot = [getAperturePhotometry(baselineImages[i], BSBX[i], BSBY[i]) for i in range(len(baselineImages))]

    BSBMean = 0
    for item in BSBPhot:
        BSBMean += item
    BSBMean /= len(BSBPhot)
    print(BSBMean)
    TSBNorm = (TSBPhot - BSBMean)/BSBMean
    plt.hist(TSBNorm, bins = 30)

    BSTMean = 0
    for item in BSTPhot:
        BSTMean += item
    BSTMean /= len(BSTPhot)
    print(BSTMean)
    TSTNorm = (TSTPhot - BSTMean)/BSTMean
    plt.hist(TSTNorm, bins = 30)
    NTSTNorm = (NTSTPhot - BSTMean)/BSTMean
    plt.hist(NTSTNorm, bins = 30)
    plt.show()

main("Images/FitsImages/Ryan1Trail1.fits","Images/NoSatellites/Ryan1Trail1_no_trail.fits","Images/NeverTrail", [3229.499999999995, 0], [2139.6999999999966, 3519], False)
