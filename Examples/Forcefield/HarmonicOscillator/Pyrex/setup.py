# Use this script to compile the Pyrex module by running
#
#         python setup.py build_ext --inplace
#
# Then run "python test.py"

from distutils.core import setup, Extension
from Pyrex.Distutils import build_ext
import os, sys

compile_args = []
include_dirs = ['../../../../Include']

from Scientific import N
try:
    num_package = N.package
except AttributeError:
    num_package = "Numeric"
if num_package == "NumPy":
    compile_args.append("-DNUMPY=1")
    import numpy.distutils.misc_util
    include_dirs.extend(numpy.distutils.misc_util.get_numpy_include_dirs())

setup (name = "HarmonicOscillatorForceField",
       version = "1.0",
       description = "Harmonic Oscillator forcefield term for MMTK",

       py_modules = ['HarmonicOscillatorFF'],
       ext_modules = [Extension('MMTK_harmonic_oscillator',
                                ['MMTK_harmonic_oscillator.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs)],
       cmdclass = {'build_ext': build_ext}
       )
