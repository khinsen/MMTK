#!/usr/bin/env python

package_name = "MMTK"
#package_name = "python-mmtk"

from distutils.core import setup, Extension
from distutils.command.sdist import sdist
from distutils.command.install_data import install_data
from distutils import dir_util
from distutils.filelist import FileList, translate_pattern
import distutils.sysconfig
sysconfig = distutils.sysconfig.get_config_vars()

import os, sys, types
from glob import glob

class Dummy:
    pass
pkginfo = Dummy()
execfile('MMTK/__pkginfo__.py', pkginfo.__dict__)


headers = glob(os.path.join ("Include", "MMTK", "*.h"))

paths = [os.path.join('MMTK', 'ForceFields'),
         os.path.join('MMTK', 'ForceFields', 'Amber'),
         os.path.join('MMTK', 'Database', 'Atoms'),
         os.path.join('MMTK', 'Database', 'Groups'),
         os.path.join('MMTK', 'Database', 'Molecules'),
         os.path.join('MMTK', 'Database', 'Complexes'),
         os.path.join('MMTK', 'Database', 'Proteins'),
         os.path.join('MMTK', 'Database', 'PDB'),
         os.path.join('MMTK', 'Tools', 'TrajectoryViewer')]
data_files = []
for dir in paths:
    files = []
    for f in glob(os.path.join(dir, '*')):
        if f[-3:] != '.py' and f[-4:-1] != '.py' and os.path.isfile(f):
            files.append(f)
    data_files.append((dir, files))


class ModifiedFileList(FileList):

    def findall(self, dir=os.curdir):
        from stat import ST_MODE, S_ISREG, S_ISDIR, S_ISLNK
        list = []
        stack = [dir]
        pop = stack.pop
        push = stack.append
        while stack:
            dir = pop()
            names = os.listdir(dir)
            for name in names:
                if dir != os.curdir:
                    fullname = os.path.join(dir, name)
                else:
                    fullname = name
                stat = os.stat(fullname)
                mode = stat[ST_MODE]
                if S_ISREG(mode):
                    list.append(fullname)
                elif S_ISDIR(mode) and not S_ISLNK(mode):
                    list.append(fullname)
                    push(fullname)
        self.allfiles = list


class modified_sdist(sdist):

    def run (self):

        self.filelist = ModifiedFileList()
        self.check_metadata()
        self.get_file_list()
        if self.manifest_only:
            return
        self.make_distribution()

    def make_release_tree (self, base_dir, files):
        self.mkpath(base_dir)
        dir_util.create_tree(base_dir, files,
                             verbose=self.verbose, dry_run=self.dry_run)
        if hasattr(os, 'link'):         # can make hard links on this system
            link = 'hard'
            msg = "making hard links in %s..." % base_dir
        else:                           # nope, have to copy
            link = None
            msg = "copying files to %s..." % base_dir
        if not files:
            self.warn("no files to distribute -- empty manifest?")
        else:
            self.announce(msg)
        for file in files:
            if os.path.isfile(file):
                dest = os.path.join(base_dir, file)
                self.copy_file(file, dest, link=link)
            elif os.path.isdir(file):
                dir_util.mkpath(os.path.join(base_dir, file))
            else:
                self.warn("'%s' not a regular file or directory -- skipping"
                          % file)

class modified_install_data(install_data):

    def run(self):
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')
        return install_data.run(self)

#################################################################
# Check various compiler/library properties

libraries = []
if sys.platform != 'win32':
    return_code = os.system(('%s config/libm_test.c -o' % sysconfig['CC']) +
                            ' config/libm_test >/dev/null 2>1')
    if return_code != 0:
        libraries.append('m')

macros = []
try:
    from Scientific.MPI import world
except ImportError:
    world = None
if world is not None:
    if type(world) == types.InstanceType:
        world = None
if world is not None:
    macros.append(('WITH_MPI', None))

if sys.platform != 'win32':
    command = '%s config/erfc_test.c ' % sysconfig['CC']
    for lib in libraries:
        command = command + '-l'+lib+' '
    command = command + '-o config/erfc_test >/dev/null 2>1'
    return_code = os.system(command)
    if return_code == 0:
        macros.append(('LIBM_HAS_ERFC', None))

command = '%s config/intlen.c -o config/intlen' % sysconfig['CC']
os.system(command)
result = os.popen('config/intlen', 'r').read().strip()
if result == '8':
    macros.append(('_LONG64_', None))

if sys.version_info[0] == 2 and sys.version_info[1] >= 2:
    macros.append(('EXTENDED_TYPES', None))

#################################################################
# System-specific optimization options

high_opt = []
if sys.platform[:5] == 'linux' and sysconfig['CC'][:3] == 'gcc':
    high_opt = ['-O3', '-ffast-math', '-fomit-frame-pointer',
                '-fkeep-inline-functions']
if sys.platform == 'darwin' and sysconfig['CC'][:3] == 'gcc':
    high_opt = ['-O3', '-ffast-math', '-fomit-frame-pointer',
                '-fkeep-inline-functions', '-falign-loops=16']
if sys.platform == 'aix4':
    high_opt = ['-O4']
if sys.platform == 'odf1V4':
    high_opt = ['-O2', '-fp_reorder', '-ansi_alias', '-ansi_args']

high_opt.append('-g')

#################################################################

setup (name = package_name,
       version = pkginfo.__version__,
       description = "Molecular Modelling Toolkit",
       long_description=
"""
The Molecular Modelling Toolkit (MMTK) is an Open Source program
library for molecular simulation applications. In addition to providing
ready-to-use implementations of standard algorithms, MMTK serves as a
code basis that can be easily extended and modified to deal with
standard and non-standard problems in molecular simulations.
""",
       author = "Konrad Hinsen",
       author_email = "hinsen@llb.saclay.cea.fr",
       url = "http://dirac.cnrs-orleans.fr/MMTK/",
       license = "LGPL",

       packages = ['MMTK', 'MMTK.ForceFields', 'MMTK.ForceFields.Amber',
                   'MMTK.NormalModes', 'MMTK.Tk', 'MMTK.Tools',
                   'MMTK.Tools.TrajectoryViewer'],
       headers = headers,
       ext_package = 'MMTK.'+sys.platform,
       ext_modules = [Extension('lapack_mmtk',
                                ['Src/lapack_subset.c', 'Src/lapack_mmtk.c'],
                                extra_compile_args = high_opt,
                                include_dirs=['Include'],
                                libraries=libraries,
                                define_macros=macros),
                      Extension('MMTK_DCD',
                                ['Src/MMTK_DCD.c', 'Src/ReadDCD.c'],
                                include_dirs=['Include'],
                                libraries=libraries,
                                define_macros=macros),
                      Extension('MMTK_deformation',
                                ['Src/MMTK_deformation.c'],
                                extra_compile_args = high_opt,
                                include_dirs=['Include'],
                                libraries=libraries,
                                define_macros=macros),
                      Extension('MMTK_dynamics',
                                ['Src/MMTK_dynamics.c'],
                                include_dirs=['Include'],
                                libraries=libraries,
                                define_macros=macros),
                      Extension('MMTK_minimization',
                                ['Src/MMTK_minimization.c'],
                                include_dirs=['Include'],
                                libraries=libraries,
                                define_macros=macros),
                      Extension('MMTK_surface',
                                ['Src/MMTK_surface.c'],
                                include_dirs=['Include'],
                                libraries=libraries,
                                define_macros=macros),
                      Extension('MMTK_trajectory',
                                ['Src/MMTK_trajectory.c'],
                                include_dirs=['Include'],
                                libraries=libraries,
                                define_macros=macros),
                      Extension('MMTK_universe',
                                ['Src/MMTK_universe.c'],
                                include_dirs=['Include'],
                                libraries=libraries,
                                define_macros=macros),
                      Extension('MMTK_forcefield',
                                ['Src/MMTK_forcefield.c',
                                 'Src/bonded.c', 'Src/nonbonded.c',
                                 'Src/ewald.c', 'Src/sparsefc.c',
                                 'Src/dpmta/mpole/mpe_fft.c',
                                 'Src/dpmta/mpole/mpe_misc.c',
                                 'Src/dpmta/mpole/mpe_mpoleC.c',
                                 'Src/dpmta/mpole/mpe_allocC.c',
                                 'Src/dpmta/mpole/mpe_mpoleLJ.c',
                                 'Src/dpmta/mpole/mpe_allocLJ.c',
                                 'Src/dpmta/src/dpmta_serial.c',
                                 'Src/dpmta/src/dpmta_slvmkcell.c',
                                 'Src/dpmta/src/dpmta_slvmcalc.c',
                                 'Src/dpmta/src/dpmta_slvpcalc.c',
                                 'Src/dpmta/src/dpmta_slvmkil.c',
                                 'Src/dpmta/src/dpmta_slvmkhl.c',
                                 'Src/dpmta/src/dpmta_slvcompute.c',
                                 'Src/dpmta/src/dpmta_slvmacro.c',
                                 'Src/dpmta/src/dpmta_slvscale.c',
                                 'Src/dpmta/src/dpmta_timer.c',
                                 'Src/dpmta/src/dpmta_slvglobals.c',
                                 'Src/dpmta/src/dpmta_distmisc.c'],
                                extra_compile_args = high_opt,
                                include_dirs=['Include', 'Src',
                                              'Src/dpmta/src',
                                              'Src/dpmta/mpole'],
                                define_macros = [('WITH_DPMTA', None),
                                                 ('SERIAL', None),
                                                 ('VIRIAL', None),
                                                 ('MACROSCOPIC', None)]
                                                + macros,
                                libraries=libraries),
                      Extension('MMTK_energy_term',
                                ['Src/MMTK_energy_term.c'],
                                include_dirs=['Include'],
                                libraries=libraries,
                                define_macros=macros),                      
                      ],

       data_files = data_files,
       scripts = ['tviewer'],

       cmdclass = {'sdist': modified_sdist,
                   'install_data': modified_install_data},
       )
