import numpy as np
import numpy.fft as F
import numpy as np

class FourierFit:
    def __init__(self, x, y, pruning):
        l = y.shape[0]
        amps = F.rfft(y)
        freqs = F.rfftfreq(l)

        self.cosFreqs = freqs
        self.sinFreqs = freqs
        
        self.cosAmps = 2*np.real(amps)/l
        self.sinAmps = -2*np.imag(amps)/l
        self.cosAmps[0] /= 2
        
        if l%2 == 0:
            self.cosAmps[-1] /= 2


        scale = (l-1)/(x[-1]-x[0])

        self.cosFreqs = self.cosFreqs*scale
        self.sinFreqs = self.sinFreqs*scale
        self.translate = -x[0]

        absCosAmps = np.abs(self.cosAmps)
        absSinAmps = np.abs(self.sinAmps)

        maskCos = (absCosAmps>0.0000001)
        maskSin = (absSinAmps>0.0000001)

        absCosAmps = absCosAmps/(self.cosFreqs+0.00001)
        absSinAmps = absSinAmps/(self.sinFreqs+0.00001)

        self.cosFreqs = self.cosFreqs[maskCos]
        self.cosAmps = self.cosAmps[maskCos]
        absCosAmps = absCosAmps[maskCos]

        self.sinFreqs = self.sinFreqs[maskSin]
        self.sinAmps = self.sinAmps[maskSin]
        absSinAmps = absSinAmps[maskSin]
        
        for i in range(0,pruning):
            if np.min(absCosAmps)<np.min(absSinAmps):
                minPos = np.argmin(absCosAmps)
                mask = np.ones(len(self.cosAmps), dtype=bool)
                mask[minPos] = False

                self.cosAmps = self.cosAmps[mask]
                absCosAmps = absCosAmps[mask]
                self.cosFreqs = self.cosFreqs[mask]
                    
            else:
                minPos = np.argmin(absSinAmps)
                mask = np.ones(len(self.sinAmps), dtype=bool)
                mask[minPos] = False
                
                self.sinAmps = self.sinAmps[mask]
                absSinAmps = absSinAmps[mask]
                self.sinFreqs = self.sinFreqs[mask]
        
        

    def __call__(self, x):
        y = 0*x
        x2 = x+self.translate
        for i in range(0,len(self.cosFreqs)):
            y += self.cosAmps[i]*np.cos(2*np.pi*self.cosFreqs[i]*x2)
            
        for i in range(0,len(self.sinFreqs)):   
            y += self.sinAmps[i]*np.sin(2*np.pi*self.sinFreqs[i]*x2)

        return y
