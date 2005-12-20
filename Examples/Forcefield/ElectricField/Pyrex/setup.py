# Use this script to compile the C module by running
#
#         python setup.py build_ext --inplace
#
# Then run "python test.py"

from distutils.core import setup, Extension
from Pyrex.Distutils import build_ext

# Patch for Pyrex.Distutils to use include paths
# This should become unnecessary in some future Pyrex release.
def pyrex_compile(self, source):
    import Pyrex.Compiler.Main
    import copy
    options = copy.deepcopy(Pyrex.Compiler.Main.default_options)
    if self.include_dirs:
        options.include_path.extend(self.include_dirs)
    result = Pyrex.Compiler.Main.compile(source, options)
    if result.num_errors <> 0:
        sys.exit(1)
build_ext.pyrex_compile = pyrex_compile
# End of Pyrex.Distutils patch

setup (name = "ElectricField",
       version = "1.0",
       description = "Electric field term for MMTK",

       py_modules = ['ElectricField'],
       ext_modules = [Extension('MMTK_electric_field',
                                ['MMTK_electric_field.pyx'])],
       cmdclass = {'build_ext': build_ext}
       )
