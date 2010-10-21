#!/usr/bin/env python

from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import numpy.distutils.misc_util
import os, sys

compile_args = ['-g', '-O0']
include_dirs = ['.', '../../Include']

from Scientific import N
assert N.package == "NumPy"

compile_args.append("-DNUMPY=1")
include_dirs.extend(numpy.distutils.misc_util.get_numpy_include_dirs())

setup (name = "MMTK-VelocityVerlet",
       version = "2.7",
       description = "Velocity Verlet for the Molecular Modelling Toolkit",
       author = "Konrad Hinsen",
       author_email = "hinsen@cnrs-orleans.fr",
       url = "http://dirac.cnrs-orleans.fr/MMTK/",
       license = "CeCILL-C",

       ext_modules = [Extension('VelocityVerlet',
                                ['./VelocityVerlet.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      Extension('VelocityVerletPI',
                                ['./VelocityVerletPI.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),],
       cmdclass = {'build_ext': build_ext},
       )
