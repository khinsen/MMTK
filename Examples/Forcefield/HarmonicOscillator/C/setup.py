# Use this script to compile the C module by running
#
#         python setup.py build_ext --inplace
#
# Then run "python test.py"

from distutils.core import setup, Extension

setup (name = "HarmonicOscillatorForceField",
       version = "1.0",
       description = "Harmonic Oscillator forcefield term for MMTK",

       py_modules = ['HarmonicOscillatorFF'],
       ext_modules = [Extension('MMTK_harmonic_oscillator',
                                ['MMTK_harmonic_oscillator.c'])]
       )
