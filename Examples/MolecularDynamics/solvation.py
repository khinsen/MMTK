# Solvation of a protein with water.
#
# The solvation procedure consists of three steps:
#
# - The universe containing the protein is scaled up, and the required
#   number of water molecules is added at random positions, but without
#   any overlap between molecules. The system is so dilute that random
#   placements are easily possible.
#
# - The universe is slowly scaled down do its original size, with
#   each scaling step followed by some energy minimization and
#   molecular dynamics steps.
#
# - A molecular dynamics run at constant pressure and temperature is
#   used to put the system into a well-defined thermodynamic state.
#

from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber94ForceField
from MMTK.Environment import NoseThermostat, AndersenBarostat
from MMTK.Trajectory import Trajectory, TrajectoryOutput, LogOutput
from MMTK.Dynamics import VelocityVerletIntegrator, VelocityScaler, \
                          BarostatReset, TranslationRemover
import MMTK.Solvation

# Create the solute.
protein = Protein('bala1')

# Put the solvent in a standard configuration: center of mass at the
# coordinate origin, principal axes of inertia parallel to the coordinate axes.
protein.normalizeConfiguration()

# Define density, pressure,  and temperature of the solvent.
water_density = 1.*Units.g/Units.cm**3
temperature = 300.*Units.K
pressure = 1.*Units.atm

# Calculate the box size as the boundary box of the protein plus an
# offset. Note: a much larger offset should be used in real applications.
box = protein.boundingBox()
box = box[1]-box[0]+Vector(0.5, 0.5, 0.5)

# Create a periodic universe. The force field is intentionally created with
# a rather small cutoff to speed up the solvation process.
universe = OrthorhombicPeriodicUniverse(tuple(box),
                                        Amber94ForceField(1., 1.))
universe.protein = protein

# Find the number of solvent molecules.
print MMTK.Solvation.numberOfSolventMolecules(universe,'water',water_density),\
      "water molecules will be added"

# Scale up the universe and add the solvent molecules.
MMTK.Solvation.addSolvent(universe, 'water', water_density)
print "Solvent molecules have been added, now shrinking universe..."

# Shrink the universe back to its original size, thereby compressing
# the solvent to its real density.
MMTK.Solvation.shrinkUniverse(universe, temperature, 'solvation.nc')
print "Universe has been compressed, now equilibrating..."

# Set a better force field and add thermostat and barostat.
#
# Note: For efficiency, optimized Ewald parameters should be used
# in a real application. The barostat relaxation time must be
# adjusted to the system size; it should be chosen smaller than
# for a realistic simulation in order to reach the final
# volume faster.
universe.setForceField(Amber94ForceField(1.4, {'method': 'ewald'}))
universe.thermostat = NoseThermostat(temperature)
universe.barostat = AndersenBarostat(pressure, 0.1*Units.ps)

# Create an integrator and a trajectory.
integrator = VelocityVerletIntegrator(universe, delta_t=1.*Units.fs)

trajectory = Trajectory(universe, "equilibration.nc", "w",
                        "Equilibration (NPT ensemble)")

# Start an NPT integration with periodic rescaling of velocities
# and resetting of the barostat. The number of steps required to
# reach a stable volume depends strongly on the system!
output_actions = [TrajectoryOutput(trajectory,
                                   ('configuration', 'energy', 'thermodynamic',
                                    'time', 'auxiliary'), 0, None, 100),
                  LogOutput("equilibration.log", ('time', 'energy'),
                            0, None, 100)]
integrator(steps = 1000,
           actions = [TranslationRemover(0, None, 200),
                      BarostatReset(0, None, 20),
                      VelocityScaler(temperature, 0., 0, None, 20)]
           + output_actions)

# Close the equilibration trajectory
trajectory.close()
