# Example: Create a cubic water box
#
# This example requires only trivial modification to perform solvatation
# of another molecule. Simply put the molecule into the universe before
# the water molecule is added, and put its bounding sphere on the list
# of excluded regions. And don't forget to calculate the final box
# size (real_size) correctly.
#

from MMTK import *
from MMTK.ForceFields import Amber94ForceField
from MMTK.Trajectory import Trajectory, TrajectoryOutput
from MMTK.Minimization import SteepestDescentMinimizer
from MMTK.Dynamics import VelocityVerletIntegrator, VelocityScaler, \
                          TranslationRemover
from MMTK.Random import randomPointInBox

# Set the number of molecules and the temperature
n_molecules = 20
temperature = 300.*Units.K

# Calculate the size of the box
density = 1.*Units.g/(Units.cm)**3
number_density = density/Molecule('water').mass()
real_size = (n_molecules/number_density)**(1./3.)

# Create the universe with a larger initial size
current_size = 5.*real_size
world = CubicPeriodicUniverse(current_size,
                              Amber94ForceField(1.2*Units.nm,
                                                {'method': 'ewald'}))

# Add solvent molecules at random positions, avoiding the neighbourhood
# of previously placed molecules
excluded_regions = []
for i in range(n_molecules):
    m = Molecule('water', position = randomPointInBox(current_size))
    while 1:
	s = m.boundingSphere()
	collision = 0
	for region in excluded_regions:
	    if s.intersectWith(region) is not None:
		collision = 1
		break
	if not collision:
	    break
	m.translateTo(randomPointInBox(current_size))
    world.addObject(m)
    excluded_regions.append(s)

# Reduce energy
minimizer = SteepestDescentMinimizer(world, step_size = 0.05*Units.Ang)
minimizer(steps = 100)

# Set velocities and equilibrate for a while
world.initializeVelocitiesToTemperature(temperature)
integrator = VelocityVerletIntegrator(world,
                                      actions=[VelocityScaler(300., 10.),
                                               TranslationRemover()])
integrator(steps = 200)

# Scale down the system in small steps
while current_size > real_size:

    scale_factor = max(0.95, real_size/current_size)
    for object in world:
	object.translateTo(scale_factor*object.position())
    current_size = scale_factor*current_size
    world.setSize(current_size)

    print 'Current size: ', current_size
    stdout.flush()

    minimizer(steps = 100)
    integrator(steps = 200)

    save(world, 'water'+`n_molecules`+'.intermediate.setup')
    
# Final equilibration
trajectory = Trajectory(world, 'water.nc', 'w', 'Final equilibration')
integrator(steps = 1000,
           actions = [TrajectoryOutput(trajectory, ("time", "energy",
                                                    "thermodynamic",
                                                    "configuration"),
                                       0, None, 10)])
trajectory.close()

# Save final system
save(world, 'water'+`n_molecules`+'.setup')
