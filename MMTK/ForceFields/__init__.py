# Force field initialization
#
# Written by Konrad Hinsen
#

"""
Force fields

:undocumented: BondedInteractions
:undocumented: ForceField
:undocumented: NonBondedInteractions
:undocumented: MMForceField
"""

__docformat__ = 'restructuredtext'

import os, string, sys

class _ForceFieldLoader(object):

    def __init__(self, module, object):
        self.module = module
        self.object = object
        self.globals = globals()

    __safe_for_unpickling__ = True

    def __call__(self, *args, **kw):
        ffc = getattr(__import__(self.module, self.globals), self.object)
        ffc.description = _description
        return apply(ffc, args, kw)

def _description(self):
    return 'ForceFields.' + self.__class__.__name__ + `self.arguments`

ff_list = open(os.path.join(os.path.split(__file__)[0],
                            'force_fields')).readlines()
for line in ff_list:
    line = string.split(line)
    exec line[0] + "=_ForceFieldLoader(line[1], line[2])"

try:
    default_energy_threads = string.atoi(os.environ['MMTK_ENERGY_THREADS'])
except KeyError:
    default_energy_threads = 1

del os
del string
del sys
del ff_list
del line

