#!/usr/bin/env python

package_name = "MMTK-LangevinDynamics"
version_number = "1.0"

from distutils.core import setup, Extension
import sys

setup (name = "MMTK-LangevinDynamics",
       version = "1.2",
       description = "Langevin dynamics module for the Molecular Modelling Toolkit",
       author = "Konrad Hinsen",
       author_email = "hinsen@llb.saclay.cea.fr",
       url = "http://dirac.cnrs-orleans.fr/MMTK/",
       license = "LGPL",

       ext_modules = [Extension('MMTK_langevin',
                                ['./MMTK_langevin.c', 'ranf.c', 'pmath_rng.c'],
                                include_dirs=['./']),
                      ],
       )
