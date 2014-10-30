# Constant-temperature constant-pressure MD simulation of Argon.
#

from MMTK import *
from MMTK.ForceFields import LennardJonesForceField
from MMTK.Environment import NoseThermostat, AndersenBarostat
from MMTK.Trajectory import Trajectory, TrajectoryOutput, LogOutput
from MMTK.Dynamics import VelocityVerletIntegrator, VelocityScaler, \
                          TranslationRemover, BarostatReset
from MMTK.ProgressOutput import ProgressOutput
import string
from Scientific.IO.TextFile import TextFile

# Open the configuration file and read the box size.
conf_file = TextFile('argon.conf.gz')
lx, ly, lz = map(string.atof, string.split(conf_file.readline()))

# Construct a periodic universe using a Lennard-Jones (noble gas) force field
# with a cutoff of 15 Angstrom.
universe = OrthorhombicPeriodicUniverse((lx*Units.Ang,
                                         ly*Units.Ang, lz*Units.Ang),
                                        LennardJonesForceField(15.*Units.Ang))
# Read the atom positions and construct the atoms.
while 1:
    line = conf_file.readline()
    if not line: break
    x, y, z = map(string.atof, string.split(line))
    universe.addObject(Atom('Ar', position=Vector(x*Units.Ang,
                                                  y*Units.Ang, z*Units.Ang)))

# Define thermodynamic parameters.
temperature = 94.4*Units.K
pressure = 1.*Units.atm

# Add thermostat and barostat.
universe.thermostat = NoseThermostat(temperature)
universe.barostat = AndersenBarostat(pressure)

# Initialize velocities.
universe.initializeVelocitiesToTemperature(temperature)

# Create trajectory and integrator.
trajectory = Trajectory(universe, "argon_npt.nc", "w", "Argon NPT test")
integrator = VelocityVerletIntegrator(universe, delta_t=10*Units.fs)

# Periodical actions for trajectory output and text log output.
output_actions = [TrajectoryOutput(trajectory,
                                   ('configuration', 'energy', 'thermodynamic',
                                    'time', 'auxiliary'), 0, None, 20),
                  LogOutput("argon.log", ('time', 'energy'), 0, None, 100)]

# Do some equilibration steps, rescaling velocities and resetting the
# barostat in regular intervals.
integrator(steps = 2000,
           actions = [TranslationRemover(0, None, 100),
                      VelocityScaler(temperature, 0.1*temperature,
                                     0, None, 100),
                      BarostatReset(100),
                      ProgressOutput(2000)] + output_actions)

# Do some "production" steps.
integrator(steps = 2000,
           actions = [TranslationRemover(0, None, 100),
                      ProgressOutput(2000)] + output_actions)

# Close the trajectory.
trajectory.close()
