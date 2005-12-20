# Force field initialization
#
# Written by Konrad Hinsen
# last revision: 2002-2-26
#

class _ForceFieldLoader:

    def __init__(self, module, object):
        self.module = module
        self.object = object
        self.globals = globals()

    __safe_for_unpickling__ = 1

    def __call__(self, *args, **kw):
        ffc = getattr(__import__(self.module, self.globals), self.object)
        return apply(ffc, args, kw)

import os, string, sys

ff_list = open(os.path.join(os.path.split(__file__)[0],
                            'force_fields')).readlines()
for line in ff_list:
    line = string.split(line)
    exec line[0] + "=_ForceFieldLoader(line[1], line[2])"

if sys.modules.has_key('pythondoc'):

    from Amber.AmberForceField import Amber94ForceField
    Amber94ForceField.__module__ = 'MMTK.ForceFields'
    from LennardJonesFF import LennardJonesForceField
    LennardJonesForceField.__module__ = 'MMTK.ForceFields'
    from DeformationFF import DeformationForceField
    DeformationForceField.__module__ = 'MMTK.ForceFields'
    from CalphaFF import CalphaForceField
    CalphaForceField.__module__ = 'MMTK.ForceFields'

try:
    default_energy_threads = string.atoi(os.environ['MMTK_ENERGY_THREADS'])
except KeyError:
    default_energy_threads = 1

del os
del string
del sys
del ff_list
del line

