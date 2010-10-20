from MMTK import *
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber99ForceField
from VelocityVerlet import VelocityVerletIntegrator
from MMTK.Dynamics import Heater, TranslationRemover, RotationRemover
from MMTK.Trajectory import Trajectory, TrajectoryOutput, \
                            RestartTrajectoryOutput, StandardLogOutput
import time

# Define system
universe = InfiniteUniverse(Amber99ForceField(mod_files=['frcmod.ff99SB']))
universe.protein = Protein('bala1')

# Initialize velocities
universe.initializeVelocitiesToTemperature(50.*Units.K)
print 'Temperature: ', universe.temperature()
print 'Momentum: ', universe.momentum()
print 'Angular momentum: ', universe.angularMomentum()

# Create integrator
integrator = VelocityVerletIntegrator(universe, delta_t=1.*Units.fs)

# Heating and equilibration
integrator(steps=1000,
                    # Heat from 50 K to 300 K applying a temperature
                    # change of 0.5 K/fs; scale velocities at every step.
	   actions=[Heater(50.*Units.K, 300.*Units.K, 0.5*Units.K/Units.fs,
                           0, None, 1),
                    # Remove global translation every 50 steps.
		    TranslationRemover(0, None, 50),
                    # Remove global rotation every 50 steps.
		    RotationRemover(0, None, 50),
                    # Log output to screen every 100 steps.
                    StandardLogOutput(100)])

# "Production" run
trajectory = Trajectory(universe, "bala1.nc", "w", "A simple test case")
integrator(steps=100,
                      # Remove global translation every 50 steps.
           actions = [TranslationRemover(0, None, 50),
                      # Remove global rotation every 50 steps.
                      RotationRemover(0, None, 50),
                      # Write every second step to the trajectory file.
                      TrajectoryOutput(trajectory, ("time", "energy",
                                                    "thermodynamic",
                                                    "configuration"),
                                       0, None, 2),
                      # Write restart data every fifth step.
                      RestartTrajectoryOutput("restart.nc", 5),
                      # Log output to screen every 10 steps.
                      StandardLogOutput(10)])
trajectory.close()

print "Starting background thread..."
thread = integrator(steps = 10000, background = True,
                    actions = [TranslationRemover(0, None, 50),
                               RotationRemover(0, None, 50)])
print "Monitoring progress:"
while thread.is_alive():
    state = thread.copyState()
    if state is not None:
        print "CPU time:", time.clock()
        print "Simulation time:", state['time']
        print "Potential energy:", state['potential_energy']
    time.sleep(0.1)
