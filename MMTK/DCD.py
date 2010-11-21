# This module implements a DCD reader/writer
#
# Written by Lutz Ehrlich
# Adapted to MMTK conventions by Konrad Hinsen


"""
Reading and writing of DCD trajectory files

The DCD format for trajectories is used by CHARMM, X-Plor,
and NAMD. It can be read by various visualization programs.

The DCD format is defined as a binary (unformatted) Fortran
format and is therefore platform-dependent.

:undocumented: writePDB
:undocumented: writeDCD
"""

__docformat__ = 'restructuredtext'

import MMTK_DCD
from MMTK import PDB, Trajectory, Units
from Scientific import N


class DCDReader(Trajectory.TrajectoryGenerator):

    """
    Reader for DCD trajectories (CHARMM/X-Plor)

    A DCDReader reads a DCD trajectory and "plays back" the
    data as if it were generated directly by an integrator.
    The universe for which the DCD file is read must be
    perfectly compatible with the data in the file, including
    an identical internal atom numbering. This can be guaranteed
    only if the universe was created from a PDB file that is
    compatible with the DCD file without leaving out any
    part of the system.

    Reading is started by calling the reader object.
    The following data categories and variables are available for
    output:

    * category "time": time
    * category "configuration": configuration
    """

    default_options = {}

    available_data = ['configuration', 'time']

    restart_data = None

    def __init__(self, universe, **options):
        """
        :param universe: the universe for which the information from the
                         trajectory file is read
        :param options: keyword options
        :keyword dcd_file: the name of the DCD trajecory file to be read
        :keyword actions: a list of actions to be executed periodically
                          (default is none)
        """
        Trajectory.TrajectoryGenerator.__init__(self, universe, options)

    def __call__(self, **options):
        self.setCallOptions(options)
        configuration = self.universe.configuration()
        MMTK_DCD.readDCD(self.universe, configuration.array,
                         self.getActions(), self.getOption('dcd_file'))

        
def writeDCD(vector_list, dcd_file_name, factor, atom_order=None,
             delta_t=0.1, conf_flag=1):
    universe = vector_list[0].universe
    npoints = universe.numberOfPoints()
    if atom_order is None:
        atom_order = N.arrayrange(npoints)
    else:
        atom_order = N.array(atom_order)

    i_start = 0       # always start at frame 0
    n_savc  = 1       # save every frame
    fd = MMTK_DCD.writeOpenDCD(dcd_file_name, npoints, len(vector_list),
                               i_start, n_savc, delta_t)
    for vector in vector_list:
        if conf_flag:
            vector = universe.contiguousObjectConfiguration(None, vector)
        array = factor*vector.array
        x = N.take(array[:, 0], atom_order).astype(N.Float16)
        y = N.take(array[:, 1], atom_order).astype(N.Float16)
        z = N.take(array[:, 2], atom_order).astype(N.Float16)
        MMTK_DCD.writeDCDStep(fd, x, y, z)
    MMTK_DCD.writeCloseDCD(fd)

def writePDB(universe, configuration, pdb_file_name):
    offset = None
    if universe is not None:
        configuration = universe.contiguousObjectConfiguration(None,
                                                               configuration)
    pdb = PDB.PDBOutputFile(pdb_file_name, 'xplor')
    pdb.write(universe, configuration)
    sequence = pdb.atom_sequence
    pdb.close()
    return sequence

def writeDCDPDB(conf_list, dcd_file_name, pdb_file_name, delta_t=0.1):
    """
    Write a sequence of configurations to a DCD file and generate
    a compatible PDB file.

    :param conf_list: the sequence of configurations
    :type conf_list: sequence of :class:`~MMTK.ParticleProperties.Configuration`
    :param dcd_file_name: the name of the DCD file
    :type dcd_file_name: str
    :param pdb_file_name: the name of the PDB file
    :type pdb_file_name: str
    :param delta_t: the time step between two configurations
    :type delta_t: float
    """
    universe = conf_list[0].universe
    sequence = writePDB(universe, conf_list[0], pdb_file_name)
    # For now, write the first bead of each atom:
    indices = [a.index for a in sequence]
    # For writing all beads, this should be:
    # indices = sum((range(a.index, a.index+a.nbeads) for a in sequence), [])
    # but this requires a change in PDB output.
    writeDCD(conf_list, dcd_file_name, 1./Units.Ang, indices, delta_t, 1)


def writeVelocityDCDPDB(vel_list, dcd_file_name, pdb_file_name, delta_t=0.1):
    """
    Write a sequence of velocity particle vectors to a DCD file and generate
    a compatible PDB file.

    :param vel_list: the sequence of velocity particle vectors
    :type vel_list: sequence of :class:`~MMTK.ParticleProperties.ParticleVector`
    :param dcd_file_name: the name of the DCD file
    :type dcd_file_name: str
    :param pdb_file_name: the name of the PDB file
    :type pdb_file_name: str
    :param delta_t: the time step between two velocity sets
    :type delta_t: float
    """
    universe = vel_list[0].universe
    sequence = writePDB(universe, universe.configuration(), pdb_file_name)
    # For now, write the first bead of each atom:
    indices = [a.index for a in sequence]
    # For writing all beads, this should be:
    # indices = sum((range(a.index, a.index+a.nbeads) for a in sequence), [])
    # but this requires a change in PDB output.
    writeDCD(vel_list, dcd_file_name, 1./(Units.Ang/Units.akma_time),
             indices, delta_t, 0)
