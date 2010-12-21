.. _Examples:

.. |C_alpha| replace:: C\ :sub:`α`

Code Examples
#############

One of the best ways to learn how to use a new tool is to look at
examples. The examples given in this manual were adapted from
real-life MMTK applications. They are also contained in the
MMTK distribution (directory "Examples") for direct use and modification.

The example molecules, system sizes, parameters, etc.,
were chosen to reduce execution time as much as possible, in order to
enable you to run the examples interactively step by step to see how
they work. If you plan to modify an example program for your own use,
don't forget to check all parameters carefully to make sure that you
obtain reasonable results.


.. _Example-MolecularDynamics: 

- Molecular Dynamics examples

  - The file
    :doc:`argon.py <Examples/MolecularDynamics/argon.py>`
    contains a simulation of liquid argon at constant temperature and
    pressure.
  - The file
    :doc:`protein.py <Examples/MolecularDynamics/protein.py>`
    contains a simulation of a small (very small) protein in vacuum.
  - The file
    :doc:`restart.py <Examples/MolecularDynamics/restart.py>`
    shows how the simulation started in
    :doc:`protein.py <Examples/MolecularDynamics/protein.py>`
    can be continued.
  - The file
    :doc:`solvation.py <Examples/MolecularDynamics/solvation.py>`
    contains the solvation of a protein by water molecules.

.. _Example-MonteCarlo:

- Monte-Carlo examples

  - The program
    :doc:`backbone.py <../Examples/MonteCarlo/backbone.py>`
    generates an ensemble of backbone configuration (C-alpha atoms only)
    for a protein.

.. _Example-Trajectories:

- Trajectory examples

  - The file
    :doc:`snapshot.py <../Examples/Trajectories/snapshot.py>`
    shows how a trajectory can be built up step by step from arbitrary
    data.
  - The file
    :doc:`dcd_import.py <../Examples/Trajectories/dcd_import.py>`
    converts a trajectory in DCD format (used by the programs CHARMM,
    X-Plor, and NAMD) to MMTK's format.
  - The file
    :doc:`dcd_export.py <../Examples/Trajectories/dcd_export.py>`
    converts an MMTK trajectory to DCD format (used by the programs CHARMM,
    X-Plor, and NAMD).
  - The file
    :doc:`trajectory_average.py <../Examples/Trajectories/trajectory_average.py>`
    calculates an average structure from a trajectory.
  - The file
    :doc:`trajectory_extraction.py <../Examples/Trajectories/trajectory_extraction.py>`
    reads a trajectory and writes a new one containing only a subset of the
    original universe.
  - The file
    :doc:`view_trajectory.py <../Examples/Trajectories/view_trajectory.py>`
    shows an animation of a trajectory, provided that an external molecule
    viewer with animation is available.
  - The file
    :doc:`calpha_trajectory.py <../Examples/Trajectories/calpha_trajectory.py>`
    shows how a much smaller |C_alpha|-only trajectory can be extracted from
    a trajectory containing one or more proteins.

.. _Example-NormalModes:

-  Normal mode examples

  - The file
    :doc:`modes.py <../Examples/NormalModes/modes.py>`
    contains a standard normal mode calculation for a small protein.
  - The file
    :doc:`constrained_modes.py <../Examples/NormalModes/constrained_modes.py>`
    contains a normal mode calculation for a small protein using a model
    in which each amino acid residue is rigid.
  - The file
    :doc:`calpha_modes.py <../Examples/NormalModes/calpha_modes.py>`
    contains a normal mode calculation for a mid-size protein using a
    |C_alpha| model and an elastic network model.
  - The file
    :doc:`harmonic_force_field.py <../Examples/NormalModes/harmonic_force_field.py>`
    contains a normal mode calculation for a protein using a detailed
    but still simple harmonic force field.

.. _Example-Proteins:

- Protein examples

  - The file
    :doc:`construction.py <../Examples/Proteins/construction.py>`
    shows some more complex examples of protein construction from PDB files.
  - The file
    :doc:`analysis.py <../Examples/Proteins/analysis.py>`
    demonstrates a few analysis techniques for comparing protein
    conformations.

.. _Example-DNA:

- DNA examples

  - The file
    :doc:`construction.py <../Examples/DNA/construction.py>`
    contains the construction of a DNA strand with a ligand.

.. _Example-Forcefield:

- Forcefield examples

  - Electric field term

    - A pure Python implementation (rather slow in general, but
      tolerable for a simple term like this one) is given in
      :doc:`Python/ElectricField.py <../Examples/Forcefield/ElectricField/Python/ElectricField.py>`.

    - A more efficient implementation has the evaluation code written
      in Cython (:doc:`Cython/MMTK_electric_field.pyx
      <../Examples/Forcefield/ElectricField/Cython/MMTK_electric_field.pyx>`)
      while the bookkeeping part remains in Python (:doc:`Cython/ElectricField.py <../Examples/Forcefield/ElectricField/Cython/ElectricField.py>`).

  - Harmonic oscillator term

    - A pure Python implementation (rather slow in general, but
      tolerable for a simple term like this one) is given in
      :doc:`Python/HarmonicOscillatorFF.py <../Examples/Forcefield/HarmonicOscillator/Python/HarmonicOscillatorFF.py>`.

    - A more efficient implementation has the evaluation code written
      in Cython (:doc:`Cython/MMTK_harmonic_oscillator.pyx
      <../Examples/Forcefield/HarmonicOscillator/Cython/MMTK_harmonic_oscillator.pyx>`)
      while the bookkeeping part remains in Python (:doc:`Cython/HarmonicOscillatorFF.py <../Examples/Forcefield/HarmonicOscillator/Cython/HarmonicOscillatorFF.py>`).

.. _Example-MPI:

- MPI examples (parallelization)

  - The file :doc:`md.py <../Examples/MPI/md.py>`
    contains a parallelized version of :doc:`solvation.py <../Examples/MolecularDynamics/solvation.py>`.

.. _Example-MDIntegrator:

- Molecular Dynamics integrators

  - The file :doc:`md.py <../Examples/MDIntegrator/VelocityVerlet.pyx>`
    illustrates how Molecular Dynamics integrators can be implemented
    in Cython.

.. _Example-LangevinDynamics:

- Langevin dynamics integrator

  - The files 
    :doc:`LangevinDynamics.py <../Examples/LangevinDynamics/LangevinDynamics.py>`
    and :doc:`MMTK_langevinmodule.c <../Examples/LangevinDynamics/MMTK_langevin.c>`
    implement a simple integrator for Langevin dynamics. It is meant as an 
    example of how to write integrators etc. in C, 
    but of course it can also be used directly.

.. _Example-Visualization:

-  Visualization examples

  - The file 
    :doc:`additional_objects.py <../Examples/Visualization/additional_objects.py>`
    describes the addition of custom graphics objects to the representation
    of a molecular system.

.. _Example-Miscellaneous:

- Micellaneous examples

  - The example
    :doc:`charge_fit.py <../Examples/Miscellaneous/charge_fit.py>`
    demonstrates fitting point charges to an electrostatic potential
    energy surface.
  - The file
    :doc:`construct_from_pdb.py <../Examples/Miscellaneous/construct_from_pdb.py>`
    shows how a universe can be built from a PDB file in such a way that
    the internal atom ordering is compatible. This is important for exchanging
    data with other programs.
  - The file
    :doc:`lattice.py <../Examples/Miscellaneous/lattice.py>`
    constructs molecules placed on a lattice.
  - The file
    :doc:`vector_field.py <../Examples/Miscellaneous/vector_field.py>`
    shows how vector fields can be used in the analysis and visualization
    of collective motions.
