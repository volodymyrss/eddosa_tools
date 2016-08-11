#! /usr/bin/env python

# System imports
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# lut2model extension module
_lut2model = Extension("_lut2model",
                   ["lut2model.i","lut2model.cxx"],
                   include_dirs = [numpy_include],
                   language = "c++",
                   )

# lut2model setup
setup(  name        = "lut2 modelling",
        description = "lut2 modelling",
        author      = "Volodymyr Savchenko",
        version     = "1.0",
        ext_modules = [_lut2model]
        )
