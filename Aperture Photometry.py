from photutils import CircularAperture as CAp
from photutils import CircularAnnulus as CAn
from photutils.aperture import aperture_photometry as ap
from photutils.utils import calc_total_error as cte

import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np


fitsStars1 = fits.open("RyanImages/Ryan1/Ryan1Trail2.fits")
data1 = fitsStars1[0].data
error1 = cte(data1, 0*data1, 76.0)

fitsStars2 = fits.open("RyanImages/Ryan1/Ryan1NoTrail2.fits")#RyanImages/Ryan1/Ryan1NoTrail2.fits")
data2 = fitsStars2[0].data
error2 = cte(data2, 0*data2, 76.0)

fitsStars4 = fits.open("RyanImages/Ryan1/Ryan1NoTrail2FWS.fits")
data4 = fitsStars4[0].data
error4 = cte(data4, 0*data4, 76.0)

fitsStars5 = fits.open("RyanImages/Ryan1/Ryan1NoTrail2NB.fits")
data5 = fitsStars5[0].data
error5 = cte(data5, 0*data5, 76.0)

fitsStars6 = fits.open("RyanImages/Ryan1/Ryan1NoTrail2NBFWS.fits")
data6 = fitsStars6[0].data
error6 = cte(data6, 0*data6, 76.0)

fitsStars7 = fits.open("RyanImages/Ryan1/Ryan1NoTrail2TS5.fits")
data7 = fitsStars7[0].data
error7 = cte(data7, 0*data7, 76.0)

fitsStars8 = fits.open("RyanImages/Ryan1/Ryan1NoTrail2TS20.fits")
data8 = fitsStars8[0].data
error8 = cte(data8, 0*data8, 76.0)

fitsStars9 = fits.open("Images/NoSatellites/L_2019-04-08_22-59-51_c_no_trail.fits")#RyanImages/Ryan1/Ryan1NoTrail2.fits")
data9 = fitsStars9[0].data
error9 = cte(data9, 0*data9, 76.0)

fitsStars3 = fits.open("Images/FitsImages/NoSatellite/L_2019-04-08_22-28-44_c.fits")
data3 = fitsStars3[0].data
error3 = cte(data3, 0*data3, 76.0)


#plt.imshow(data1)
#plt.show()
#plt.imshow(data2)
#plt.show()
#plt.imshow(data3)
#plt.show()

stars = [(4186,1933), (4134,1921), (3970,1902), (3809,1888), (3758,1871), (2990,1772), (2343,1681), (1406,1571), (899,1512)]
stars2 = [(item[0]-41,item[1]-18) for item in stars]
aperture = CAp(stars, r=8)
aperture2 = CAp(stars2, r=8)
annulus1 = CAn(stars, r_in = 12, r_out = 18)
annulus2 = CAn(stars2, r_in = 12, r_out = 18)


aperture_masks1 = aperture.to_mask()
aperture_masks2 = aperture2.to_mask()
star_apertures1 = [aperture_masks1[i].multiply(data1) for i in range(0,len(aperture_masks1))]
star_apertures2 = [aperture_masks1[i].multiply(data2) for i in range(0,len(aperture_masks1))]
star_apertures3 = [aperture_masks2[i].multiply(data3) for i in range(0,len(aperture_masks2))]



#annulus_masks = annulus.to_mask(method = 'center')
#star_annulus = [annulus_masks[i].multiply(data2) for i in range(0,len(annulus_masks))]

#for i in range(0,len(star_apertures1)):
#    fig, axes  = plt.subplots(1,3)
#    axes[0].imshow(star_apertures1[i])
#    axes[1].imshow(star_apertures2[i])
#    axes[2].imshow(star_apertures3[i])
#    plt.show()

phot_table1 = ap(data1, [aperture,annulus1], error = error1)
phot_table2 = ap(data2, [aperture,annulus1], error = error2)
phot_table4 = ap(data4, [aperture,annulus1], error = error4)
phot_table5 = ap(data5, [aperture,annulus1], error = error5)
phot_table6 = ap(data6, [aperture,annulus1], error = error6)
phot_table7 = ap(data7, [aperture,annulus1], error = error7)
phot_table8 = ap(data8, [aperture,annulus1], error = error8)
phot_table9 = ap(data9, [aperture,annulus1], error = error9)
phot_table3 = ap(data3, [aperture2,annulus2], error = error3)

print(phot_table1)

phot1 = phot_table1["aperture_sum_0"]/aperture.area#-aperture.area*phot_table1["aperture_sum_1"]/annulus1.area
phot2 = phot_table2["aperture_sum_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus1.area
phot4 = phot_table4["aperture_sum_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area
phot5 = phot_table5["aperture_sum_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area
phot6 = phot_table6["aperture_sum_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area
phot7 = phot_table7["aperture_sum_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area
phot8 = phot_table8["aperture_sum_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area
phot9 = phot_table9["aperture_sum_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area
phot3 = phot_table3["aperture_sum_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area

phot1error = phot_table1["aperture_sum_err_0"]/aperture.area#-aperture.area*phot_table1["aperture_sum_1"]/annulus1.area
phot2error = phot_table2["aperture_sum_err_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus1.area
phot4error = phot_table4["aperture_sum_err_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area
phot5error = phot_table5["aperture_sum_err_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area
phot6error = phot_table6["aperture_sum_err_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area
phot7error = phot_table7["aperture_sum_err_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area
phot8error = phot_table8["aperture_sum_err_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area
phot9error = phot_table9["aperture_sum_err_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area
phot3error = phot_table3["aperture_sum_err_0"]/aperture.area#-aperture2.area*phot_table1["aperture_sum_1"]/annulus2.area

print(phot1)
print(phot2)
print(phot3)

phot1Diff = 100*(phot1-phot3)#/phot3
phot2Diff = 100*np.abs(phot2-phot3)#/phot3
phot4Diff = np.abs(phot4-phot3)/phot3
phot5Diff = np.abs(phot5-phot3)/phot3
phot6Diff = np.abs(phot6-phot3)/phot3
phot7Diff = np.abs(phot7-phot3)/phot3
phot8Diff = np.abs(phot8-phot3)/phot3
phot9Diff = 100*(phot9-phot3)#/phot3

diff1Error = 100*np.sqrt(phot1**2+phot3**2)
diff2Error = 100*np.sqrt(phot2**2+phot3**2)
diff9Error = 100*np.sqrt(phot9**2+phot3**2)

print("Average error with trail: "+str(np.mean(phot1Diff)))
print("Average error without trail: "+str(np.mean(phot2Diff)))
print("Average error without 0.5 trail: "+str(np.mean(phot9Diff)))
print("Error std with trail: " + str(np.std(phot1Diff)))
print("Error std without trail: " + str(np.std(phot2Diff)))
print("Error std without 0.5 trail: " + str(np.std(phot9Diff)))

plt.plot(np.abs(phot1Diff),np.abs(phot9Diff),"ro")
plt.plot([0,max(max(phot1Diff),max(phot9Diff))],[0,max(max(phot1Diff),max(phot9Diff))])
plt.show()

fig, axes = plt.subplots(1,2)
axes[0].errorbar([0,1,2,3,4,5,6,7,8],np.abs(np.array(phot1Diff)), yerr = np.array(diff1Error), label = "Error with trail")
axes[0].errorbar([0,1,2,3,4,5,6,7,8],np.abs(np.array(phot9Diff)), yerr = np.array(diff9Error), label = "Error without trail")
axes[0].plot(list(np.abs(phot1Diff)), "b:")
axes[0].plot(list(np.abs(phot9Diff)), "g--")
axes[0].set_title("Photometry Comparison With Trail And Trail Removed")
axes[0].set_xlabel("Star Number")
axes[0].set_ylabel("Photometry Error")
axes[0].legend()

axes[1].errorbar([0,1,2,3,4,5,6,7,8],np.array(phot1Diff), yerr = np.array(diff1Error), label = "Error with trail")
axes[1].errorbar([0,1,2,3,4,5,6,7,8],np.array(phot9Diff), label = "Error without trail")
axes[1].plot(list(phot1Diff), "b:")
axes[1].plot(list(phot9Diff), "g--")
axes[1].set_title("Photometry Comparison With Trail And Trail Removed")
axes[1].set_xlabel("Star Number")
axes[1].set_ylabel("Photometry Error")
axes[1].legend()
plt.show()

