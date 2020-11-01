from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    #ext_modules = cythonize("Main.py")
    #    ext_modules = cythonize("TrailRemovalCode/PiecewiseCubic.pyx"), include_dirs=[numpy.get_include()]
    ext_modules = cythonize("TrailRemovalCode/FourierFit.pyx"), include_dirs=[numpy.get_include()]
)
