import numpy as np
from numpy import load
import matplotlib.pyplot as plt
import TrailRemovalCode.FourierFit as FF

pars = load('pars.npy')
offset= load('offset.npy')

print(4550/pars.shape[0])

f = FF.FourierFit(pars,offset,0)
plt.plot(pars,offset)
plt.plot(pars,f(pars))
plt.show()

for i in range(4500,4700,10):
    print(i)
    f = FF.FourierFit(pars,offset,i)
    plt.plot(pars,offset)
    plt.plot(pars,f(pars))
    plt.show()
