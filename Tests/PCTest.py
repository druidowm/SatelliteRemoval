from TrailRemovalCode.PiecewiseCubic import PiecewiseCubic as PC
import numpy as np
from matplotlib import pyplot as plt

x = np.arange(0.0,10)
y = np.arange(0.0,10)
x,y = np.meshgrid(x,y)
z = x**2+y**2

pc = PC(z)

x1 = np.arange(0.0,10)
y1 = np.arange(0.0,10)
x1,y1 = np.meshgrid(x1,y1)
z1 = x1*0

for i in range(0,len(x1)):
    for j in range(0,len(x1[0])):
        z1[i][j] = pc(x1[i][j],y1[i][j])


fig,axes = plt.subplots(2,1)
axes[0].imshow(z)
axes[1].imshow(z1)
plt.show()
