# Use this script to compile the C module by running
#
#         python setup.py build_ext --inplace
#
# Then run "python test.py"

from distutils.core import setup, Extension

setup (name = "ElectricField",
       version = "1.0",
       description = "Electric field term for MMTK",

       py_modules = ['ElectricField'],
       ext_modules = [Extension('MMTK_electric_field',
                                ['MMTK_electric_field.c'])]
       )
