import numpy as np
import math

class PiecewiseCubic:
    def __init__(self, img):
        x = np.arange(0,4,1)
        y = np.arange(0,4,1)
        x,y = np.meshgrid(x,y)
        x = x.ravel()
        y = y.ravel()
        zeros = np.ones(x.shape[0])
        x2 = x**2
        y2 = y**2
        x3 = x**3
        y3 = y**3
        xy = x*y
        x2y = x2*y
        xy2 = y2*x
        x3y = x3*y
        xy3 = x*y3
        x2y2 = x2*y2
        x2y3 = x2*y3
        x3y2 = x3*y2
        x3y3 = x3*y3
        A = np.array([zeros,
                                  x, y, x2, y2,
                                  x3, y3, xy,
                                  x2y, xy2,
                                  x2y2, x3y,
                                  xy3, x3y2,
                                  x2y3, x3y3]).T

        self.functions = [[None for j in range(0,img.shape[1]-3)] for i in range(0,img.shape[0]-3)]
        for i in range(0,img.shape[0]-3):
            for j in range(0,img.shape[1]-3):
                mask = img[i:i+4,j:j+4]

                mask = mask.ravel()

                coeff = np.linalg.solve(A, mask)
                self.functions[i][j] = coeff

        self.shape = img.shape



    def __call__(self,x,y):
        if x<2:
            i = 0
        elif 2<=x<self.shape[1]-2:
            i = int(x-2)
        else:
            i = self.shape[1]-4

        if y<2:
            j = 0
        elif 2<=y<self.shape[0]-2:
            j = int(y-2)
        else:
            j = self.shape[0]-4

        coeffs = self.functions[j][i]

        xNew = x-i
        yNew = y-j
        x2 = xNew**2
        y2 = yNew**2
        x3 = xNew**3
        y3 = yNew**3

        val = coeffs[0]
        val += coeffs[1]*xNew
        val += coeffs[2]*yNew
        val += coeffs[3]*x2
        val += coeffs[4]*y2
        val += coeffs[5]*(xNew**3)
        val += coeffs[6]*(yNew**3)
        val += coeffs[7]*(xNew*yNew)
        val += coeffs[8]*x2*yNew
        val += coeffs[9]*y2*xNew
        val += coeffs[10]*x2*y2
        val += coeffs[11]*x3*yNew
        val += coeffs[12]*y3*xNew
        val += coeffs[13]*x3*y2
        val += coeffs[14]*x2*y3
        val += coeffs[15]*x3*y3

        return val

        
