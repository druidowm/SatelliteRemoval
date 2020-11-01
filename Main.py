import TrailRemovalCode.RemoveTrail as RT
from TrailRemovalCode.PiecewiseCubic import PiecewiseCubic as PC
import StarRemovalCode.StarRemoval as SR
import Solve
import numpy as np
import astropy as ap
from astropy.io import fits
from matplotlib import cm
from matplotlib import pyplot as plt

def solve():
    dPath = "Images/FitsImages"
    outPath = "Images/NoSatellites"
    Solve.batchSolveWithClient(dPath, outPath, 10, showProgress = True, showNoStars = True, showTrail = True, showNoTrail = True, showNoTrailWithTrail = True, inverted = True)

def testSolve():
    dPath = "Images/FitsImages/Ryan1Trail1.fits"
    outPath = "Images/NoSatellites/NOTRAIL.fit"
    #Solve.solveTest(dPath, outPath, 10, [0, 1402.599999999994], [4655, 1991.299999999997])
    Solve.solveTest(dPath, outPath, 10, [3229.499999999995, 0], [2139.6999999999966, 3519], inverted = True)
    #Solve.solveTest(dPath, outPath, 10, [2173.6999999999985, 0], [4655, 2270.9999999999936])
    #Solve.solveTest(dPath, outPath, 10, [45,2001], [3013, 186], inverted = True)


#imgPath = "Images/FitsImages/Starlink2.fit"
#Solve.solveWithClient(imgPath, 10, showProgress = True, showNoStars = False, showTrail = True, showNoTrail = True, showNoTrailWithTrail = False, inverted = False)
solve()
#testSolve()