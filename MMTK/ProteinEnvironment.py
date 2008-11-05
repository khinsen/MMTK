# This module defines the environment in which protein definition files
# are executed.

from MMTK.Proteins import PeptideChain
from MMTK.Collections import Collection
from MMTK.PDB import PDBConfiguration
from MMTK.Units import *
from Scientific.Geometry import Vector
from Scientific.Geometry.Transformation import Translation, Rotation
from copy import copy, deepcopy
