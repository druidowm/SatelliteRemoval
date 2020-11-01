import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt

class CubicFit:
    def __init__(self, img):
        x = np.arange(0,img.shape[1],1)
        y = np.arange(0,img.shape[0],1)
        x,y = np.meshgrid(x,y)
        x = x.ravel()
        y = y.ravel()
        imgF = img.ravel()
        zeros = np.ones(x.shape[0])
        x2 = x**2
        y2 = y**2
        x3 = x**3
        y3 = y**3
        xy = x*y
        x2y = x2*y
        xy2 = y2*x
        #x3y = x3*y
        #xy3 = x*y3
        x2y2 = x2*y2
        #x2y3 = x2*y3
        #x3y2 = x3*y2
        #x3y3 = x3*y3
        A = np.array([zeros,
                        x, y, x2, y2,#x3, y3,
                        xy,
                        x2y, xy2,
                        x2y2,# x3y,
                        #xy3, x3y2,x2y3, x3y3
                          ]).T
        
        self.coeffs,a,b,c = np.linalg.lstsq(A,imgF)

    def __call__(self,x,y):
        x2 = x**2
        y2 = y**2
        #x3 = x**3
        #y3 = y**3

        val = self.coeffs[0]
        val += self.coeffs[1]*x
        val += self.coeffs[2]*y
        val += self.coeffs[3]*x2
        val += self.coeffs[4]*y2
        #val += self.coeffs[5]*(x**3)
        #val += self.coeffs[6]*(y**3)
        val += self.coeffs[5]*(x*y)
        val += self.coeffs[6]*x2*y
        val += self.coeffs[7]*y2*x
        val += self.coeffs[8]*x2*y2
        #val += self.coeffs[11]*x3*y
        #val += self.coeffs[12]*y3*x
        #val += self.coeffs[13]*x3*y2
        #val += self.coeffs[14]*x2*y3
        #val += self.coeffs[15]*x3*y3

        return val
