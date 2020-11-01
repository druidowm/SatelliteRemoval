from photutils import CircularAperture as CAp
from photutils import CircularAnnulus as CAn
from photutils.aperture import aperture_photometry as ap
from photutils.utils import calc_total_error as cte

import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

class Photometry:
    def __init__(self, imgTrail, imgNoTrail, compImg, gain):
        fitsTrail = fits.open(imgTrail)
        self.dataTrail = fitsTrail[0].data
        self.errorTrail = cte(self.dataTrail, 0*self.dataTrail, gain)

        fitsNoTrail = fits.open(imgNoTrail)
        self.dataNoTrail = fitsNoTrail[0].data
        self.errorNoTrail = cte(self.dataNoTrail, 0*self.dataNoTrail, gain)

        fitsComp = fits.open(compImg)
        self.dataComp = fitsComp[0].data
        self.errorComp = cte(self.dataComp, 0*self.dataComp, gain)
        print(self.errorComp)



    def getApertures(self, stars1, offset, rAp, annulus = False, rAnIn = None, rAnOut = None, show = False):
        stars2 = [(item[0]+offset[0],item[1]+offset[1]) for item in stars1]

        self.annulus = annulus
        
        self.apertureTrail = CAp(stars1, r=rAp)
        self.apertureComp = CAp(stars2, r=rAp)

        if annulus:
            self.annulusTrail = CAn(stars1, r_in = rAnIn, r_out = rAnOut)
            self.annulusComp = CAn(stars2, r_in = rAnIn, r_out = rAnOut)

        if show:
            apertureMaskTrail = self.apertureTrail.to_mask(method = 'center')
            starTrail = [apertureMaskTrail[i].multiply(self.dataTrail) for i in range(0,len(apertureMaskTrail))]
            
            starNoTrail = [apertureMaskTrail[i].multiply(self.dataNoTrail) for i in range(0,len(apertureMaskTrail))]

            apertureMaskComp = self.apertureComp.to_mask(method = 'center')
            starComp = [apertureMaskComp[i].multiply(self.dataComp) for i in range(0,len(apertureMaskComp))]
            
            if annulus:
                annulusMaskTrail = self.annulusTrail.to_mask(method = 'center')
                backTrail = [annulusMaskTrail[i].multiply(self.dataTrail) for i in range(0,len(annulusMaskTrail))]
                
                backNoTrail = [annulusMaskTrail[i].multiply(self.dataNoTrail) for i in range(0,len(annulusMaskTrail))]

                annulusMaskComp = self.annulusComp.to_mask(method = 'center')
                backComp = [annulusMaskComp[i].multiply(self.dataComp) for i in range(0,len(annulusMaskComp))]

                for i in range(0,len(starTrail)):
                    fig, axes  = plt.subplots(2,3)
                    
                    axes[0][0].imshow(starTrail[i])
                    axes[0][0].set_title("Star with Trail")
                    
                    axes[0][1].imshow(starNoTrail[i])
                    axes[0][1].set_title("Star with Trail Removed")
                    
                    axes[0][2].imshow(starComp[i])
                    axes[0][2].set_title("Star without Trail")
                    
                    axes[1][0].imshow(backTrail[i])
                    axes[1][0].set_title("Background with Trail")
                    
                    axes[1][1].imshow(backNoTrail[i])
                    axes[1][1].set_title("Background with Trail Removed")
                    
                    axes[1][2].imshow(backComp[i])
                    axes[1][2].set_title("Background without Trail")
                    
                    plt.show()


            else:
                for i in range(0,len(starTrail)):
                    fig, axes  = plt.subplots(1,3)
                    
                    axes[0].imshow(starTrail[i])
                    axes[0].set_title("Star with Trail")
                    
                    axes[1].imshow(starNoTrail[i])
                    axes[1].set_title("Star with Trail Removed")
                    
                    axes[2].imshow(starComp[i])
                    axes[2].set_title("Star without Trail")
                    
                    plt.show()



    def getPhotometry(self):
        if self.annulus:
            trail = ap(self.dataTrail, [self.apertureTrail,self.annulusTrail], error = self.errorTrail)
            noTrail = ap(self.dataNoTrail, [self.apertureTrail,self.annulusTrail], error = self.errorNoTrail)
            comp = ap(self.dataComp, [self.apertureComp,self.annulus1Comp], error = self.errorComp)

            self.photTrail = trail["aperture_sum_0"]/self.apertureTrail.area - trail["aperture_sum_1"]/self.annulusTrail.area
            self.photNoTrail = noTrail["aperture_sum_0"]/self.apertureTrail.area - noTrail["aperture_sum_1"]/self.annulusTrail.area
            self.photComp = comp["aperture_sum_0"]/self.apertureComp.area - comp["aperture_sum_1"]/self.annulusComp.area

            self.photTrailError = np.sqrt((trail["aperture_sum_err_0"]/self.apertureTrail.area)**2+(trail["annulus_sum_err_0"]/self.annulusTrail.area)**2)
            self.photNoTrailError = np.sqrt((noTrail["aperture_sum_err_0"]/self.apertureTrail.area)**2+(noTrail["annulus_sum_err_0"]/self.annulusTrail.area)**2)
            self.photCompError = np.sqrt((comp["aperture_sum_err_0"]/self.apertureComp.area)**2+(comp["annulus_sum_err_0"]/self.annulusComp.area)**2)

        else:
            trail = ap(self.dataTrail, self.apertureTrail, error = self.errorTrail)
            noTrail = ap(self.dataNoTrail, self.apertureTrail, error = self.errorNoTrail)
            comp = ap(self.dataComp, self.apertureComp, error = self.errorComp)

            self.photTrail = trail["aperture_sum"]/self.apertureTrail.area
            self.photNoTrail = noTrail["aperture_sum"]/self.apertureTrail.area
            self.photComp = comp["aperture_sum"]/self.apertureComp.area

            self.photTrailError = trail["aperture_sum_err"]/self.apertureTrail.area
            self.photNoTrailError = noTrail["aperture_sum_err"]/self.apertureTrail.area
            self.photCompError = comp["aperture_sum_err"]/self.apertureComp.area

        print(self.photCompError)




    def getPhotometryDiffs(self, ab = False, printOut = False, show = False):
        if ab:
            self.photDiffsTrail = np.abs(self.photTrail - self.photComp)
            self.photDiffsNoTrail = np.abs(self.photNoTrail - self.photComp)

        else:
            self.photDiffsTrail = self.photTrail - self.photComp
            self.photDiffsNoTrail = self.photNoTrail - self.photComp


        self.photDiffsErrorTrail = np.sqrt(self.photTrailError**2+self.photCompError**2)
        self.photDiffsErrorNoTrail = np.sqrt(self.photNoTrailError**2+self.photCompError**2)
        print(self.photDiffsErrorTrail)


    def getPhotometryPercentDiffs(self):
        self.photPercentDiffsTrail = self.photDiffsTrail/self.photComp
        self.photPercentDiffsNoTrail = self.photDiffsNoTrail/self.photComp

        self.photPercentDiffsErrorTrail = self.photPercentDiffsTrail*np.sqrt((self.photDiffsErrorTrail/self.photDiffsTrail)**2 + (self.photCompError/self.photComp)**2)
        self.photPercentDiffsErrorNoTrail = self.photPercentDiffsNoTrail*np.sqrt((self.photDiffsErrorNoTrail/self.photDiffsNoTrail)**2 + (self.photCompError/self.photComp)**2)


    def printInfo(self):
        print("Average error with trail: "+str(np.mean(self.photDiffsTrail)))
        print("Average error without trail: "+str(np.mean(self.photDiffsNoTrail)))
        print("Error std with trail: " + str(np.std(self.photDiffsTrail)))
        print("Error std without trail: " + str(np.std(self.photDiffsNoTrail)))


    def plotError(self):
        x = [i for i in range(0,len(self.photDiffsTrail))]
        
        plt.errorbar(x,np.array(self.photDiffsTrail), yerr = np.array(self.photDiffsErrorTrail), label = "Error with trail")
        plt.errorbar(x,np.array(self.photDiffsNoTrail), yerr = np.array(self.photDiffsErrorNoTrail), label = "Error without trail")
        #axes[0].plot(list(np.abs(phot1Diff)), "b:")
        #axes[0].plot(list(np.abs(phot9Diff)), "g--")
        plt.title("Photometry Comparison With Trail And Trail Removed")
        plt.xlabel("Star Number")
        plt.ylabel("Photometry Error")
        plt.legend()
        plt.show()


#p = Photometry("RyanImages/Ryan1/Ryan1Trail2.fits","Images/NoSatellites/L_2019-04-08_22-59-51_c_no_trail.fits","Images/FitsImages/NoSatellite/L_2019-04-08_22-28-44_c.fits", 76.0)
p = Photometry("RyanImages/Ryan1/Ryan1Trail2.fits","RyanImages/Ryan1/Ryan1NoTrail2.fits","Images/FitsImages/NoSatellite/L_2019-04-08_22-28-44_c.fits", 76.0)
p.getApertures([(4186,1933), (4134,1921), (3970,1902), (3809,1888), (3758,1871), (2990,1772), (2343,1681), (1406,1571), (899,1512)], [-41, -18], 8, show = True)
p.getPhotometry()
p.getPhotometryDiffs()
p.plotError()
