# This program demonstrates how point charges can be fitted to
# known values of the electrostatic potential energy surface.
# In a real application, these values would come from an
# ab-initio calculation.
#
from MMTK import *
from MMTK.ChargeFit import ChargeFit, TotalChargeConstraint, \
                           evaluationPoints

# Construct a small system of two atoms.
a1 = Atom('C', position=Vector(-0.05,0.,0.))
a2 = Atom('C', position=Vector( 0.05,0.,0.))
system = Collection(a1, a2)

# Define charges. Note that their sum is not zero.
a1.charge = -0.75
a2.charge = 0.76

# Calculate the electrostatic potential energy of the two charges.
# This is obviously unrealistic; in a real application, the potential
# energy surface would not be due to known point charges, but come
# from ab-initio calculations.
points = []
for r in evaluationPoints(system, 50):
    p = 0.
    for atom in system.atomList():
        p = p + atom.charge/(r-atom.position()).length()
    points.append((r, p*Units.electrostatic_energy))

# Define a charge constraint that ensures neutrality of the total system.
constraints = [TotalChargeConstraint(system, 0.)]

# Perform the fit and print the two fitted charges. They are not equal
# to the input charges due to the neutrality constraint.
fit = ChargeFit(system, points, constraints)
print fit[a1], fit[a2]
