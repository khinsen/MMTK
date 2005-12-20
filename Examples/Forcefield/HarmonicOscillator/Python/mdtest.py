from MMTK import *
from HarmonicOscillatorFF import HarmonicOscillatorForceField
from MMTK.Dynamics import VelocityVerletIntegrator
from MMTK.Trajectory import TrajectoryOutput

universe = InfiniteUniverse()
universe.atom1 = Atom('C', position=Vector(0., 0., 1.05))
universe.atom2 = Atom('C', position=Vector(0., 1.05, 0.))

ff1 = HarmonicOscillatorForceField(universe.atom1, Vector(0., 0., 1.), 100.)
ff2 = HarmonicOscillatorForceField(universe.atom2, Vector(0., 1., 0.), 100.)
universe.setForceField(ff1+ff2)

universe.initializeVelocitiesToTemperature(300.*Units.K)
integrator = VelocityVerletIntegrator(universe, delta_t = 1.0*Units.fs)
integrator(steps = 1000)
