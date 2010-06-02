#!/usr/bin/env python

from distutils.core import setup, Extension
import os, sys

compile_args = []
include_dirs = ['.']

from Scientific import N
try:
    num_package = N.package
except AttributeError:
    num_package = "Numeric"
if num_package == "NumPy":
    compile_args.append("-DNUMPY=1")
    if sys.platform == 'win32':
        include_dirs.append(os.path.join(sys.prefix,
                            "Lib/site-packages/numpy/core/include"))
    else:
        include_dirs.append(os.path.join(sys.prefix,
                            "lib/python%s.%s/site-packages/numpy/core/include"
                             % sys.version_info [:2]))

setup (name = "MMTK-LangevinDynamics",
       version = "2.7",
       description = "Langevin dynamics module for the Molecular Modelling Toolkit",
       author = "Konrad Hinsen",
       author_email = "hinsen@cnrs-orleans.fr",
       url = "http://dirac.cnrs-orleans.fr/MMTK/",
       license = "CeCILL-C",

       py_modules = ['LangevinDynamics'],
       ext_modules = [Extension('MMTK_langevin',
                                ['./MMTK_langevin.c', 'ranf.c', 'pmath_rng.c'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      ],
       )
