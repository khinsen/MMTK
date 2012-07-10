MMTK User's Guide
#################
for MMTK |version|

by Konrad Hinsen

.. |C_alpha| replace:: C\ :sub:`α`

Introduction
############

The Molecular Modelling Toolkit (MMTK)
presents a new approach to molecular simulations. It is not a
"simulation program" with a certain set of functions that can be used
by writing more or less flexible "input files", but a collection of
library modules written in an easy-to-learn high-level programming
language, `Python <http://www.python.org>`_. This approach
offers three important advantages: 

- Application programs can use the full power of a general and well-designed
  programming language.

- Application programs can profit from the large set of other
  libraries that are available for Python. These may be scientific or
  non-scientific; for example, it is very easy to write simulation or
  analysis programs with graphical user interfaces (using the
  module Tkinter in the Python standard library), or couple scientific
  calculations with a Web server.

- Any user can provide useful additions in separate modules, whereas adding
  features to a monolithic program requires at least the cooperation of
  the original author.

To further encourage collaborative code development, MMTK uses a very
unrestrictive licensing policy. MMTK is free software, just like
Python. Although MMTK is copyrighted, anyone is allowed to use it for
any purpose, including commercial ones, as well as to modify and
redistribute it (a more precise description is given in the copyright
statement that comes with the code).

This manual describes version |version| of MMTK. The 2.x versions contain
some incompatible changes with respect to earlier versions (1.x), most
importantly a package structure that reduces the risk of name
conflicts with other Python packages, and facilitates future
enhancements. There are also many new features and improvements to
existing functions.

Using MMTK requires a basic knowledge of object-oriented programming
and Python. Newcomers to this subject should have a look at the
introductory section in this manual and at the `Python tutorial
<http://www.python.org/doc/tut/tut.html>`_ (which also comes with the
Python interpreter). There are also numerous `books on Python
<http://www.python.org/doc/Books.html>`_ that are useful in getting
started. Even without MMTK, Python is a very useful programming
language for scientific use, allowing rapid development and testing
and easy interfacing to code written in low-level languages such as
Fortran or C.

This manual consists of several introductory chapters and a
:ref:`module reference <Reference>`. The introductory chapters explain
how common tasks are handled with MMTK, but they do not describe all
of its features, nor do they contain a full documentation of functions
or classes. This information can be found in the
:ref:`module -reference <Reference>`, which describes all classes and functions
intended for end-user applications module by module, using
documentation extracted directly from the source code.  References
from the introductory sections to the module reference facilitate
finding the relevant documentation.

Overview
########

This chapter explains the basic structure of MMTK and its view of
molecular systems. Every MMTK user should read it at least once.

Using MMTK
==========

MMTK applications are ordinary Python programs, and can be written
using any standard text editor. For interactive use it is recommended
to use one of the many tools available for interactive Python
programming.

MMTK tries to be as user-friendly as possible for interactive use. For
example, lengthy calculations can usually be interrupted by typing
Control-C. This will result in an error message ("Keyboard
Interrupt"), but you can simply go on typing other commands.
Interruption is particularly useful for energy minimization and
molecular dynamics: you can interrupt the calculation at any time,
look at the current state or do some analysis, and then continue.

Modules
=======

MMTK is a package consisting of various modules, most of them written
in Python, and some in C for efficiency. The individual modules are
described in the :ref:`module -reference <Reference>`. The basic
definitions that almost every application needs are collected in the
top-level module, MMTK. The first line of most applications is
therefore::

    from MMTK import *

The definitions that are specific to particular applications reside in
submodules within the package MMTK. For example, force fields are
defined in :mod:`MMTK.ForceFields`, and peptide
chain and protein objects are defined in :mod:`MMTK.Proteins`.

Python provides two ways to access objects in modules and submodules.
The first one is importing a module and referring to objects in it,
e.g.::

    import MMTK
    import MMTK.ForceFields
    universe = MMTK.InfiniteUniverse(MMTK.ForceFields.Amber99ForceField())

The second method is importing all or some objects *from*
a module::

    from MMTK import InfiniteUniverse
    from MMTK.ForceFields import Amber99ForceField
    universe = InfiniteUniverse(Amber99ForceField())

These two import styles can also be mixed according to convience.
In order to prevent any confusion, all objects are referred to by
their full names in this manual. The Amber force field object
is thus called :class:`MMTK.ForceFields.Amber99ForceField`.
Of course the user is free to use selective imports in order to
be able to use such objects with shorter names.

Objects
=======

MMTK is an object-oriented system.  Since objects are everywhere and
everything is an object, it is useful to know the most important
object types and what can be done with them. All object types in MMTK
have meaningful names, so it is easy to identify them in practice. The
following overview contains only those objects that a user will see
directly. There are many more object types used by MMTK internally,
and also some less common user objects that are not mentioned here.

Chemical objects
----------------

These are the objects that represent the parts of a molecular system:

- atoms
- groups
- molecules
- molecular complexes

These objects form a simple hierarchy: complexes consist of
molecules, molecules consist of groups and atoms, groups consist of
smaller groups and atoms. All of these, except for groups,
can be used directly to construct a molecular system. Groups can
only be used in the definitions of other groups and molecules in the
:ref:`overview-database`.

A number of operations can be performed on chemical objects, which
can roughly be classified into inquiry (constituent atoms, bonds, center
of mass etc.) and modification (translate, rotate).

There are also specialized versions of some of these objects. For example,
MMTK defines proteins as special complexes, consisting of peptide
chains, which are special molecules. They offer a range of special
operations (such as selecting residues or constructing the positions of
missing hydrogen atoms) that do not make sense for molecules in
general.

Collections
-----------

Collection objects represent arbitrary collections of chemical
objects.  They are used to be able to refer to collections as single
entities. For example, you might want to call all water molecules
collectively "solvent". Most of the operations on chemical objects
are also available for collections.

Force fields
------------

Force field objects represent a precise description of force
fields, i.e. a complete recipe for calculating the potential energy
(and its derivatives) for a given molecular system. In other words,
they specify not only the functional form of the various interactions,
but also all parameters and the prescriptions for applying these
parameters to an actual molecular system.

Universes
---------

Universes define complete molecular systems, i.e. they contain
chemical objects. In addition, they describe interactions within the
system (by a force field), boundary conditions, external fields,
etc. Many of the operations that can be used on chemical objects can
also be applied to complete universes.

Minimizers and integrators
--------------------------

A minimizer object is a special "machine" that can find local minima
in the potential energy surface of a universe. You may consider this a
function, if you wish, but of course functions are just special
objects. Similarly, an integrator is a special "machine" that can
determine a dynamical trajectory for a system on a given potential
energy surface.

Trajectories
------------

Minimizers and integrators can produce trajectories, which are special
files containing a sequence of configurations and/or other related
information. Of course trajectory objects can also be read for
analysis.

Variables
---------

Variable objects (not to be confused with standard Python variables)
describe quantities that have a value for each atom in a system, for
example positions, masses, or energy gradients. Their most common use
is for storing various configurations of a system.

Normal modes
------------

Normal mode objects contain normal mode frequencies and atomic
displacements for a given universe. MMTK provides three kinds of
normal modes which correspond to three different physical situations:

- Energetic modes represent the principal axes of the potential
  energy surface. They each have a force constant that is smallest
  for the most collective motions and highest for the most localized
  motions. Energetic modes are appropriate for describing the
  potential energy surface without reference to any particular
  dynamics, e.g. in flexibility analysis or Monte-Carlo sampling.

- Vibrational modes represent the vibrational motions of a system
  that have well-defined frequencies. Their use implies that the
  system follows Newtonian dynamics or quantum dynamics. This is
  appropriate for small molecules, and for the most localized motions
  of macromolecules.

- Brownian modes represent diffusional motions of an overdamped
  system. They describe systems with Brownian dynamics and are
  appropriate for describing the slow motions of macromolecules.

Non-MMTK objects
----------------

An MMTK application program will typically also make use of objects
provided by Python or Python library modules. A particularly useful
library is the package Scientific, which
is also used by MMTK itself. The most important objects are

- numbers (integers, real number, complex numbers), provided by Python

- vectors (in 3D coordinate space) provided by the module Scientific.Geometry.

- character strings, provided by Python

- files, provided by Python

Of course MMTK applications can make use of the Python standard
library or any other Python modules. For example, it is possible
to write a simulation program that provides status reports via an
integrated Web server, using the Python standard module SimpleHTTPServer.

.. _overview-database:

The chemical database
=====================

For defining the chemical objects described above, MMTK uses a
database of descriptions. There is a database for atoms, one for
groups, etc. When you ask MMTK to make a specific chemical object, for
example a water molecule, MMTK looks for the definition of water in
the molecule database. A database entry contains everything there
is to know about the object it defines: its constituents and their
names, configurations, other names used e.g. for I/O, and all
information force fields might need about the objects.

MMTK comes with database entries for many common objects (water,
amino acids, etc.). For other objects you will have to write the definitions
yourself. as described in the section :ref:`database`.

Force fields
============

MMTK contains everything necessary to use the 
`Amber 99 force field <http://ambermd.org/#ff>`_ 
on proteins, DNA, and
water molecules. It uses the standard Amber parameter and modification
file format. In addition to the Amber force field, there is a simple
Lennard-Jones force field for noble gases, and a deformation force
field for normal mode calculations on large proteins.

MMTK was designed to make the addition of force field terms and the
implementation of other force fields as easy as possible. Force field
terms can be defined in Python (for ease of implementation) or in
Cython, C, or Fortran (for efficiency). This is described in the
developer's guide.

Units
=====

Since MMTK is not a black-box program, but a modular library,
it is essential for it to use a consistent unit system in which, for
example, the inverse of a frequency is a time, and the product of
a mass and the square of a velocity is an energy, without additional
conversion factors. Black-box programs can (and usually do) use
a consistent unit system internally and convert to "conventional"
units for input and output.

The unit system of MMTK consists mostly of SI units of appropriate
magnitude for molecular systems:

=============== =========
**Measurement** **Units**
=============== =========
Length          nm
--------------- ---------
Time            ps
--------------- ---------
Mass            amu (g/mol)
--------------- ---------
Energy          kJ/mol
--------------- ---------
Frequency       THz (1/ps)
--------------- ---------
Temperature     K
--------------- ---------
Charge          e
=============== =========

The module :mod:`MMTK.Units` contains convenient
conversion constants for the units commonly used in computational
chemistry. For example, a length of 2 Ångström can be
written as ``2*Units.Ang``, and a frequency can be
printed in wavenumbers with ``print frequency/Units.invcm``.

A simple example
================

The following simple example shows how a typical MMTK application
might look like. It constructs a system consisting of a single water
molecule and runs a short molecular dynamics trajectory. There are
many alternative ways to do this; this particular one was chosen
because it makes each step explicit and clear. The individual steps
are explained in the remaining chapters of the manual.
::

    # Import the necessary MMTK definitions.
    from MMTK import *
    from MMTK.ForceFields import Amber99ForceField
    from MMTK.Trajectory import Trajectory, TrajectoryOutput, StandardLogOutput
    from MMTK.Dynamics import VelocityVerletIntegrator
    # Create an infinite universe (i.e. no boundaries, non-periodic).
    universe = InfiniteUniverse(Amber99ForceField())
    # Create a water molecule in the universe.
    # Water is defined in the database.
    universe.molecule = Molecule('water')
    # Generate random velocities.
    universe.initializeVelocitiesToTemperature(300*Units.K)
    # Create an integrator.
    integrator = VelocityVerletIntegrator(universe) 
    # Generate a trajectory 
    trajectory = Trajectory(universe, "water.nc", "w") 
    # Run the integrator for 50 steps of 1 fs, printing time and energy  
    # every fifth step and writing time, energy, temperature, and the positions  
    # of all atoms to the trajectory at each step.  
    t_actions = [StandardLogOutput(5), \
                 TrajectoryOutput("time", "energy", \
                 "thermodynamic", "configuration"), 0, None, 1]
    integrator(delta_t = 1.*Units.fs, steps = 50, actions = t_actions)  
    # Close the trajectory  
    trajectory.close()  

Constructing a molecular system
###############################

The construction of a complete system for simulation or analysis
involves some or all of the following operations:

- Creating molecules and other chemical objects.

- Defining the configuration of all objects.

- Defining the "surroundings" (e.g. boundary conditions).

- Choosing a force field.

MMTK offers a large range of functions to deal with these tasks.

Creating chemical objects
=========================

Chemical objects (atoms, molecules, complexes) are created from
definitions in the :ref:`database`. Since
these definitions contain most of the necessary information, the
subsequent creation of the objects is a simple procedure.

All objects are created by their class name
(:class:`MMTK.ChemicalObjects.Atom`, :class:`MMTK.ChemicalObjects.Molecule`,
and :class:`MMTK.ChemicalObjects.Complex`) with the name
of the definition file as first parameter. Additional optional parameters
can be specified to modify the object being created. The following optional
parameters can be used for all object types:

- name=string
  Specifies a name for the object. The default name is the one given in
  the definition file.

- position=vector
  Specifies the position of the center of mass. The default is the origin.

- configuration=string
  Indicates a configuration from the configuration dictionary in the
  definition file. The default is 'default' if such an entry exists in the
  configuration dictionary. Otherwise the object is created without atomic
  positions.

Some examples with additional explanations for specific types:

- ``Atom('C')`` creates a carbon atom.

- ``Molecule('water', position=Vector(0.,0.,1.))``
  creates a water molecule using configuration 'default' and moves the
  center of mass to the indicated position.

Proteins, peptide chains, and nucleotide chains
-----------------------------------------------

MMTK contains special support for working with proteins, peptide
chains, and nucleotide chains. As described in the chapter
:ref:`database`, proteins can be described by a special database
definition file. However, it is often simpler to create macromolecular
objects directly in an application program. The classes are
:class:`MMTK.Proteins.PeptideChain`, :ref:`Class:MMTK.Proteins.Protein`,
and :class:`MMTK.NucleicAcids.NucleotideChain`.

Proteins can be created from definition files in the database,
from previously constructed peptide chain objects, or directly
from PDB files if no special manipulations are necessary.

Examples: 

- Protein('insulin') creates a protein object for 
  insulin from a database file.

- Protein('1mbd.pdb') creates a protein object for myoglobin 
  directly from a PDB file, but leaving out the
  heme group, which is not a peptide chain.

Peptide chains are created from a sequence of residues, which can be 
a :class:`MMTK.PDB.PDBPeptideChain` object, a list of three-letter
residue codes, or a string containing one-letter residue codes. In the
last two cases the atomic positions are not defined. MMTK provides
several models for the residues which provide different levels of
detail: an all-atom model, a model without hydrogen atoms, two models
containing only polar hydrogens (using different definitions of polar
hydrogens), and a model containing only the |C_alpha| atoms, with each
|C_alpha| atom having the mass of the entire residue. The last model
is useful for conformational analyses in which only the backbone
conformations are important.

The construction of nucleotide chains is very similar. The residue
list can be either a :class:`MMTK.PDB.PDBNucleotideChain` object or a
list of two-letter residue names. The first letter of a residue name
indicates the sugar type ('R' for ribose and'D' for desoxyribose), and
the second letter defines the base ('A', 'C', and'G', plus 'T' for DNA
and'U' for RNA). The models are the same as for peptide chains, except
that the |C_alpha| model does not exist.

Most frequently proteins and nucleotide chains are created from a PDB
file. The PDB files often contain solvent (water) as well, and perhaps
some other molecules. MMTK provides convenient functions for extracting
information from PDB files and for building molecules from them in the
module :mod:`MMTK.PDB`. The first step is the creation of a 
:class:`MMTK.PDB.PDBConfiguration` object from the PDB file:

::

    from MMTK.PDB import PDBConfiguration
    configuration = PDBConfiguration('some_file.pdb')

The easiest way to generate MMTK objects for all molecules in the
PDB file is then
::

    molecules = configuration.createAll()

The result is a collection of molecules, peptide chains, and
nucleotide chains, depending on the contents of the PDB files.
There are also methods for modifying the PDBConfiguration before
creating MMTK objects from it, and for creating objects
selectively. See the documentation for the modules :class:`MMTK.PDB` and Scientific.IO.PDB for details,
as well as the :ref:`Proteins <Example-Proteins>` and :ref:`DNA <Example-DNA>` examples.

Lattices
--------

Sometimes it is necessary to generate objects (atoms or molecules)
positioned on a lattice. To facilitate this task, MMTK defines lattice
objects which are essentially sequence objects containing points or
objects at points. Lattices can therefore be used like lists with
indexing and for-loops. The lattice classes
are :class:`MMTK.Geometry.RhombicLattice`, 
:class:`MMTK.Geometry.BravaisLattice`, and 
:class:`MMTK.Geometry.SCLattice`.

Random numbers
--------------

The Python standard library and the Numerical Python package provide
random number generators, and more are available in seperate packages.
MMTK provides some convenience functions that return more specialized
random quantities: random points in a universe, random velocities,
random particle displacement vectors, random orientations. These
functions are defined in module :class:`MMTK.Random`.

Collections
-----------

Often it is useful to treat a collection of several objects as a
single entity.  Examples are a large number of solvent molecules
surrounding a solute, or all sidechains of a protein. MMTK has special
collection objects for this purpose, defined as
:class:`MMTK.Collections.Collection`. Most of the methods
available for molecules can also be used on collections.

A variant of a collection is the partitioned collection, implemented in
class :class:`MMTK.Collections.PartitionedCollection`. This class
acts much like a standard collection, but groups its elements by
geometrical position in small sub-boxes. As a consequence, some geometrical
algorithms (e.g. pair search within a cutoff) are much faster, but
other operations become somewhat slower.

Creating universes
------------------

A universe describes a complete molecular system consisting of any
number of chemical objects and a specification of their interactions
(i.e. a force field) and surroundings: boundary conditions, external
fields, thermostats, etc. The universe classes are defined in module
MMTK:

- :class:`MMTK.Universe.InfiniteUniverse` represents an infinite universe,
  without any boundary or periodic boundary conditions.

- :class:`MMTK.Universe.ParallelepipedicPeriodicUniverse` represents
  a general periodic universe defined by three basis vectors.

- :class:`MMTK.Universe.OrthorhombicPeriodicUniverse` represents a periodic
  universe with an orthorhombic elementary cell, whose size is defined
  by the three edge lengths. The edges are oriented along the axes of
  the coordinate system.

- :class:`MMTK.Universe.CubicPeriodicUniverse` is a special case
  of :class:`MMTK.Universe.OrthorhombicPeriodicUniverse` in which the
  elementary cell is cubic.

Universes are created empty; the contents are then added to them.
Three types of objects can be added to a universe: chemical objects
(atoms, molecules, etc.), collections, and environment objects
(thermostats etc.). It is also possible to remove objects from a
universe.

Force fields
------------

MMTK comes with several force fields, and permits the definition of
additional force fields. Force fields are defined in module
:class:`MMTK.ForceFields.ForceField`. The most import built-in force
field is the `Amber 99 force field
<http://ambermd.org/#ff>`_, represented by the
class
:class:`MMTK.ForceFields.Amber.AmberForceField.Amber99ForceField`. It
offers several strategies for electrostatic interactions, including
Ewald summation, a fast multipole method :ref:`[Rankin2002] <Rankin2002>`, and cutoff
with charge neutralization and optional screening :ref:`[Wolf1999] <Wolf1999>`..

In addition to the Amber 99 force field, there is the older Amber 94
forcefield, a Lennard-Jones force field for noble gases (Class
:class:`MMTK.ForceFields.LennardJonesFF`) and an Elastic Network Model
force field for protein normal mode calculations
(:class:`MMTK.ForceFields.DeformationFF.CalphaForceField`).

Referring to objects and parts of objects
=========================================

Most MMTK objects (in fact all except for atoms) have a hierarchical
structure of parts of which they consist. For many operations it is
necessary to access specific parts in this hierarchy.

In most cases, parts are attributes with a specific name. For example,
the oxygen atom in every water molecule is an attribute with the name
"O". Therefore if ``w`` refers to a water molecule, then ``w.O``
refers to its oxygen atom. For a more complicated example, if ``m``
refers to a molecule that has a methyl group called "M1", then
``m.M1.C`` refers to the carbon atom of that methyl group. The names
of attributes are defined in the database.

Some objects consist of parts that need not have unique names, for
example the elements of a collection, the residues in a peptide chain,
or the chains in a protein. Such parts are accessed by indices; the
objects that contain them are Python sequence types. Some examples:

- Asking for the number of items: if ``c``
  refers to a collection, then ``len(c)`` is
  the number of its elements.

- Extracting an item: if ``p`` refers to a
  protein, then ``p[0]`` is its first peptide
  chain.

- Iterating over items: if ``p`` refers to a
  peptide chain, then ``for residue in p: print
  residue.position()`` will print the center of mass positions of
  all its residues.

Peptide and nucleotide chains also allow the operation of slicing: if
``p`` refers to a peptide chain, then ``p[1:-1]`` is a subchain
extending from the second to the next-to-last residue.

The structure of peptide and nucleotide chains
----------------------------------------------

Since peptide and nucleotide chains are not constructed from an explicit
definition file in the database, it is not evident where their
hierarchical structure comes from. But it is only the top-level
structure that is treated in a special way. The constituents of peptide
and nucleotide chains, residues, are normal group objects. The
definition files for these group objects are in the MMTK standard
database and can be freely inspected and even modified or overriden by
an entry in a database that is listed earlier in MMTKDATABASE.

Peptide chains are made up of amino acid residues, each of which is a
group consisting of two other groups, one being called "peptide" and
the other "sidechain". The first group contains the peptide group and
the C and H atoms; everything else is contained in the sidechain. The
C atom of the fifth residue of peptide chain ``p`` is therefore
referred to as ``p[4].peptide.C_alpha``.

Nucleotide chains are made up of nucleotide residues, each of which is
a group consisting of two or three other groups. One group is called
"sugar" and is either a ribose or a desoxyribose group, the second one
is called "base" and is one the five standard bases. All but the first
residue in a nucleotide chain also have a subgroup called "phosphate"
describing the phosphate group that links neighbouring residues.

Analyzing and modifying atom properties
=======================================

General operations
------------------

Many inquiry and modification operations act at the atom level and can
equally well be applied to any object that is made up of atoms,
i.e. atoms, molecules, collections, universes, etc.  These operations
are defined once in a :term:`Mix-in class` called
:class:`MMTK.Collections.GroupOfAtoms`, but are available for all
objects for which they make sense. They include inquiry-type functions
(total mass, center of mass, moment of inertia, bounding box, total
kinetic energy etc.), coordinate modifications (translation, rotation,
application of :ref:`transformation`) and coordinate comparisons (RMS
difference, optimal fits).

.. _transformation:

Coordinate transformations
--------------------------

The most common coordinate manipulations involve translations and
rotations of specific parts of a system. It is often useful to refer
to such an operation by a special kind of object, which permits the
combination and analysis of transformations as well as its application
to atomic positions.

Transformation objects specify a general displacement consisting of a
rotation around the origin of the coordinate system followed by a
translation. They are defined in the module ``Scientific.Geometry``,
but for convenience the module ``MMTK`` contains a reference to them
as well.  Transformation objects corresponding to pure translations
can be created with ``Translation(displacement)``; transformation
objects describing pure rotations with ``Rotation(axis, angle)`` or
``Rotation(rotation_matrix)``.  Multiplication of transformation
objects returns a composite transformation.

The translational component of any transformation can be obtained by
calling the method ``translation()``; the rotational component is
obtained analogously with ``rotation()``.  The displacement vector for
a pure translation can be extracted with the method
``displacement()``, a tuple of axis and angle can be extracted from a
pure rotation by calling ``axisAndAngle()``.

.. _atom_property:

Atomic property objects
-----------------------

Many properties in a molecular system are defined for each individual
atom: position, velocity, mass, etc. Such properties are represented
in special objects, defined in module MMTK:
:class:`MMTK.ParticleProperties.ParticleScalar` for scalar quantities,
:class:`MMTK.ParticleProperties.ParticleVector` for vector quantities,
and :class:`MMTK.ParticleProperties.ParticleTensor` for rank-2
tensors.  All these objects can be indexed with an atom object to
retrieve or change the corresponding value. Standard arithmetic
operations are also defined, as well as some useful methods.

Configurations
--------------

A configuration object, represented by the class
:class:`MMTK.ParticleProperties.Configuration` is a special variant of
a :class:`MMTK.ParticleProperties.ParticleVector` object.  In addition
to the atomic coordinates of a universe, it stores geometric
parameters of a universe that are subject to change, e.g. the edge
lengths of the elementary cell of a periodic universe.  Every universe
has a current configuration, which is what all operations act on by
default. It is also the configuration that is updated by
minimizations, molecular dynamics, etc. The current configuration can
be obtained by calling the method ``configuration()``.

There are two ways to create configuration objects: by making a copy
of the current configuration (with ``universe.copyConfiguration()``,
or by reading a configuration from a :ref:`trajectory <Trajectories>`.

Minimization and Molecular Dynamics
###################################

.. _Trajectories:

Trajectories
============

Minimization and dynamics algorithms produce sequences of configurations
that are often stored for later analysis. In fact, they are often the
most valuable result of a lengthy simulation run. To make sure that the
use of trajectory files is not limited by machine compatibility, MMTK
stores trajectories in `netCDF <http://www.unidata.ucar.edu/packages/netcdf/>`_
files. These files contain binary data, minimizing disk space usage, but
are freely interchangeable between different machines. In addition,
there are a number of programs that can perform standard operations on
arbitrary netCDF files, and which can therefore be used directly on MMTK
trajectory files. Finally, netCDF files are self-describing, i.e.
contain all the information needed to interpret their contents.
An MMTK trajectory file can thus be inspected and processed without
requiring any further information.

For illustrations of trajectory operations, see the :ref:`examples
<Example-Trajectories>`.

Trajectory file objects are represented by the class
:class:`MMTK.Trajectory.Trajectory`. They can be opened for reading,
writing, or modification. The data in trajectory files can be stored
in single precision or double precision; single-precision is usually
sufficient, but double-precision files are required to reproduce a
given state of the system exactly.

A trajectory is closed by calling the method ``close()``.
If anything has been written to a trajectory, closing it is required to
guarantee that all data has been written to the file. Closing a
trajectory after reading is recommended in order to prevent memory
leakage, but is not strictly required.

Newly created trajectories can contain all
objects in a universe or any subset; this is useful for limiting the
amount of disk space occupied by the file by not storing uninteresting
parts of the system, e.g. the solvent surrounding a protein. It is
even possible to create a trajectory for a subset of the atoms in a
molecule, e.g. for only the |C_alpha| atoms of a protein. The universe
description that is stored in the trajectory file contains all
chemical objects of which at least one atom is represented.

When a trajectory is opened for reading, no universe object needs
to be specified. In that case, MMTK creates a universe from
the description contained in the trajectory file. This universe will
contain the same objects as the one for which the trajectory file was
created, but not necessarily have all the properties of the original
universe (the description contains only the names and types of the
objects in the universe, but not, for example, the force field). The
universe can be accessed via the attribute universe
of the trajectory.

If the trajectory was created with partial data for some of the objects,
reading data from it will set the data for the missing parts to
"undefined". Analysis operations on such systems must be done very
carefully. In most cases, the trajectory data will contain the atomic
configurations, and in that case the "defined" atoms can be extracted
with the method ``atomsWithDefinedPositions()``.

MMTK trajectory files can store various data: atomic positions,
velocities, energies, energy gradients etc. Each trajectory-producing
algorithm offers a set of quantities from which the user can choose what
to put into the trajectory. Since a detailed selection would be
tedious, the data is divided into classes, e.g. the class "energy"
stands for potential energy, kinetic energy, and whatever other
energy-related quantities an algorithm produces.

For optimizing I/O efficiency, the data layout in a trajectory file
can be modified by the block_size parameter. Small
block sizes favour reading or writing all data for one time step,
whereas large block sizes (up to the number of steps in the trajectory)
favour accessing a few values for all time steps, e.g. scalar
variables like energies or trajectories for individual atoms. The
default value of the block size is one.

Every trajectory file contains a history of its creation. The creation
of the file is logged with time and date, as well as each operation that
adds data to it with parameters and the time/date of start and end. This
information, together with the comment and the number of atoms and steps
contained in the file, can be obtained with the function 
:func:`MMTK.Trajectory.trajectoryInfo`.

It is possible to read data from a trajectory file that is being
written to by another process. For efficiency, trajectory data is not
written to the file at every time step, but only approximately every
15 minutes. Therefore the amount of data available for reading may be
somewhat less than what has been produced already.

Options for minimization and dynamics
=====================================

Minimizers and dynamics integrators accept various optional parameter
specifications. All of them are selected by keywords, have reasonable
default values, and can be specified when the minimizer or integrator
is created or when it is called. In addition to parameters that are
specific to each algorithm, there is a general parameter ``actions``
that specifies actions that are executed periodically, including
trajectory and console output.

Periodic actions
----------------

Periodic actions are specified by the keyword parameter ``actions``
whose value is a list of periodic actions, which defaults to an empty
list. Some of these actions are applicable to any
trajectory-generating algorithm, especially the output actions. Others
make sense only for specific algorithms or specific universes,
e.g. the periodic rescaling of velocities during a Molecular Dynamics
simulation.

Each action is described by an action object. The step numbers for
which an action is executed are specified by three parameters. The
parameter ``first`` indicates the number of the first step for which
the action is executed, and defaults to 0. The parameter ``last``
indicates the last step for which the action is executed, and default
to ``None``, meaning that the action is executed indefinitely. The
parameter ``skip`` speficies how many steps are skipped between two
executions of the action. The default value of 1 means that the action
is executed at each step. Of course an action object may have
additional parameters that are specific to its action.

The output actions are defined in the module :mod:`MMTK.Trajectory`
and can be used with any trajectory-generating algorithm. They are:

- :class:`MMTK.Trajectory.TrajectoryOutput` for writing data to a
  trajectory. Note that it is possible to use several trajectory
  output actions simultaneously to write to multiple trajectories. It
  is thus possible, for example, to write a short dense trajectory
  during a dynamics run for analyzing short-time dynamics, and
  simultaneously a long-time trajectory with a larger step spacing,
  for analyzing long-time dynamics.

- :class:`MMTK.Trajectory.RestartTrajectoryOutput`, which is a
  specialized version of :class:`MMTK.Trajectory.TrajectoryOutput`.
  It writes the data that the algorithm needs in order to be restarted
  to a restart trajectory file.  A restart trajectory is a trajectory
  that stores a fixed number of steps which are reused cyclically,
  such that it always contain the last few steps of a trajectory.

- :class:`MMTK.Trajectory.LogOutput` for text output
  of data to a file.

- :class:`MMTK.Trajectory.StandardLogOutput`, a specialized
  version of :class:`MMTK.Trajectory.LogOutput` that
  writes the data classes "time" and "energy" during the whole
  simulation run to standard output.

The other periodic actions are meaningful only for Molecular Dynamics
simulations:

- :class:`MMTK.Dynamics.VelocityScaler` is used for
  rescaling the velocities to force the kinetic energy to the value
  defined by some temperature. This is usually done during initial
  equilibration.

- :class:`MMTK.Dynamics.BarostatReset` resets the
  barostat coordinate to zero and is during initial equilibration
  of systems in the NPT ensemble.

- :class:`MMTK.Dynamics.Heater` rescales the velocities
  like  :class:`MMTK.Dynamics.VelocityScaler`, but
  increases the temperature step by step.

- :class:`MMTK.Dynamics.TranslationRemover` subtracts
  the global translational velocity of the system from all individual
  atomic velocities. This prevents a slow but systematic energy flow
  into the degrees of freedom of global translation, which occurs
  with most MD integrators due to non-perfect conservation of momentum.

- :class:`MMTK.Dynamics.RotationRemover` subtracts
  the global angular velocity of the system from all individual
  atomic velocities. This prevents a slow but systematic energy flow
  into the degrees of freedom of global rotation, which occurs
  with most MD integrators due to non-perfect conservation of angular
  momentum.

Fixed atoms
-----------

During the course of a minimization or molecular dynamics algorithm,
the atoms move to different positions. It is possible to exclude
specific atoms from this movement, i.e. fixing them at their initial
positions.  This has no influence whatsoever on energy or force
calculations; the only effect is that the atoms' positions never
change. Fixed atoms are specified by giving them an attribute fixed
with a value of one. Atoms that do not have an attributefixed, or one
with a value of zero, move according to the selected algorithm.

.. _energy_minimization:

Energy minimization
===================

MMTK has two energy minimizers using different algorithms: steepest
descent (:class:`MMTK.Minimization.SteepestDescentMinimizer`) and
conjugate gradient (:class:`MMTK.Minimization.ConjugateGradientMinimizer`)
. Steepest descent minimization is very inefficient if the goal is to
find a local minimum of the potential energy. However, it has the
advantage of always moving towards the minimum that is closest to the
starting point and is therefore ideal for removing bad contacts in a
unreasonably high energy configuration. For finding local minima, the
conjugate gradient algorithm should be used.

Both minimizers accept three specific optional parameters:

- ``steps`` (an integer) to specify the maximum number of
  steps (default is 100)

- ``step_size`` (a number)
  to specify an initial step length used in the search for a minimum
  (default is 2 pm)

- ``convergence`` (a number)
  to specify the gradient norm (more precisely the root-mean-square
  length) at which the minimization should stop (default is 0.01
  kJ/mol/nm)

There are three classes of trajectory data: "energy" includes the
potential energy and the norm of its gradient, "configuration" stands
for the atomic positions, and "gradients" stands for the energy
gradients at each atom position.

The following example performs 100 steps of steepest descent
minimization without producing any trajectory or printed output:
::

    from MMTK import *
    from MMTK.ForceFields import Amber99ForceField
    from MMTK.Minimization import SteepestDescentMinimizer
    universe = InfiniteUniverse(Amber99ForceField())
    universe.protein = Protein('insulin')
    minimizer = SteepestDescentMinimizer(universe)
    minimizer(steps = 100)

See also the example file :ref:`modes.py <Example-NormalModes>`.

Molecular dynamics
==================

The techniques described in this section are illustrated by several
:ref:`examples <Example-MolecularDynamics>`.

Velocities
----------

The integration of the classical equations of motion for an atomic
system requires not only positions, but also velocities for all atoms.
Usually the velocities are initialized to random values drawn from a
normal distribution with a variance corresponding to a certain
temperature. This is done by calling the method
:func:`~MMTK.Universe.Universe.initializeVelocitiesToTemperature`
on a universe. Note that the velocities are assigned atom by atom; no
attempt is made to remove global translation or rotation of the total
system or any part of the system.

During equilibration of a system, it is common to multiply all
velocities by a common factor to restore the intended temperature. This
can done explicitly by calling the method
:func:`~MMTK.Universe.Universe.scaleVelocitiesToTemperature`
on a universe, or by using the action object :class:`MMTK.Dynamics.VelocityScaler`.

Distance constraints
--------------------

A common technique to eliminate the fastest (usually uninteresting)
degrees of freedom, permitting a larger integration time step,
is the use of distance constraints on some or all chemical bonds.
MMTK allows the use of distance constraints on any pair of
atoms, even though constraining anything but chemical bonds
is not recommended due to considerable modifications of the
dynamics of the system :ref:`[vanGunsteren1982] <vanGunsteren1982>`,
:ref:`[Hinsen1995] <Hinsen1995>`.

MMTK permits the definition of distance constraints on all atom pairs
in an object that are connected by a chemical bond by calling the
method setBondConstraints. Usually this is called
for a complete universe, but it can also be called for a chemical
object or a collection of chemical objects. The 
:func:`~MMTK.ChemicalObjects.ChemicalObject.removeDistanceConstraints`
removes all distance constraints from the object for which it is called.

Constraints defined as described above are automatically taken into
account by Molecular Dynamics integrators. It is also possible to
enforce the constraints explicitly by calling the method
:func:`~MMTK.Universe.Universe.enforceConstraints` for a universe. This has the
effect of modifying the configuration and the velocities (if
velocities exist) in order to make them compatible with the
constraints.

Thermostats and barostats
-------------------------

A standard Molecular Dynamics integration allows time averages
corresponding to the NVE ensemble, in which the number of molecules,
the system volume, and the total energy are constant. This ensemble
does not represent typical experimental conditions very well.
Alternative ensembles are the NVT ensemble, in which the temperature
is kept constant by a thermostat, and the NPT ensemble, in which
temperature and pressure are kept constant by a thermostat and a
barostat. To obtain these ensembles in MMTK, thermostat and barostat
objects must be added to a universe. In the presence of these objects,
the Molecular Dynamics integrator will use the extended-systems method
for producing the correct ensemble. The classes to be used are
:class:`MMTK.Environment.NoseThermostat` and :class:`MMTK.Environment.AndersenBarostat`.

Integration
-----------

A Molecular Dynamics integrator based on the "Velocity Verlet"
algorithm :ref:`[Swope1982] <Swope1982>`, which was extended
to handle distance constraints as well as thermostats and
barostats :ref:`[Kneller1996] <Kneller1996>`, is implemented by the
class :class:`MMTK.Dynamics.VelocityVerletIntegrator`.
It has two optional keyword parameters:

- ``steps`` (an integer) to specify the number
  of steps (default is 100)

- ``delta_t`` (a number) to specify the time step
  (default 1 fs)

There are three classes of trajectory data: "energy" includes the
potential energy and the kinetic energy, as well as the energies of
thermostat and barostat coordinates if they exist, "time" stands for the time,
"thermodynamic" stand for temperature and pressure,
"configuration" stands for the atomic positions, "velocities" stands for
the atomic velocities, and "gradients" stands for the energy gradients
at each atom position.

The following example performs a 1000 step dynamics integration, storing
every 10th step in a trajectory file and removing the total translation
and rotation every 50th step:
::

    from MMTK import *
    from MMTK.ForceFields import Amber99ForceField
    from MMTK.Dynamics import VelocityVerletIntegrator, \
                              TranslationRemover, \
                              RotationRemover
    from MMTK.Trajectory import TrajectoryOutput
    universe = InfiniteUniverse(Amber99ForceField())
    universe.protein = Protein('insulin')
    universe.initializeVelocitiesToTemperature(300.*Units.K)
    actions = [TranslationRemover(0, None, 50), \
               RotationRemover(0, None, 50), \
               TrajectoryOutput("insulin.nc", \
               ("configuration", "energy", "time"), \
               0, None, 10)]
    integrator = VelocityVerletIntegrator(universe, delta_t = 1.*Units.fs, \
                                          actions = actions)
    integrator(steps = 1000)

Snapshots
=========

A snapshot generator allows writing the current system state to a
trajectory. It works much like a zero-step minimization or dynamics run,
i.e. it takes the same optional arguments for specifying the trajectory
and protocol output. A snapshot generator is created using the
class :class:`MMTK.Trajectory.SnapshotGenerator`.

Normal modes
############

Normal mode analysis provides an analytic description of the dynamics
of a system near a minimum using an harmonic approximation to the
potential. Before a normal mode analysis can be started, the system
must be brought to a local minimum of the potential energy by
:ref:`energy minimization <energy_minimization>`, except when special
force fields designed only for normal mode analysis are used
(e.g. :class:`MMTK.ForceFields.CalphaFF.CalphaForceField`). See
also the :ref:`Normal Modes <Example-NormalModes>` examples.

A standard normal mode analysis is performed by creating a normal
modes object, implemented in the classes
:class:`MMTK.NormalModes.EnergeticModes.EnergeticModes`,
:class:`MMTK.NormalModes.VibrationalModes.VibrationalModes`, and
:class:`MMTK.NormalModes.BrownianModes.BrownianModes`. A normal
mode object behaves like a sequence of mode objects, which store the
atomic displacement vectors corresponding to each mode and the
associated force constant, vibrational frequency, or inverse
relaxation time.

For short-ranged potentials, it is advantageous to store the second
derivatives of the potential in a sparse-matrix form and to use
an iterative method to determine some or all modes. This permits
the treatment of larger systems that would normally require huge
amounts of memory. 

Another approach to deal with large systems is the restriction to
low-frequency modes which are supposed to be well representable by
linear combinations of a given set of basis vectors. The basis vectors
can be obtained from a basis for the full Cartesian space by
elimination of known fast degrees of freedom (e.g. bonds); the
module :mod:`MMTK.Subspace` contains support classes for this
approach. It is also possible to construct a suitable basis vector set
from small-deformation vector fields
(e.g. :class:`MMTK.FourierBasis.FourierBasis`).

Analysis operations
###################

Analysis is the most non-standard part of molecular simulations.
The quantities that must be calculated depend strongly on the
system and the problem under study. MMTK provides a wide range
of elementary operations that inquire the state of the system,
as well as several more complex analysis tools. Some of them are
demonstrated in the :ref:`examples` section.

Properties of chemical objects and universes
============================================

Many operations access and modify various properties of an object. They
are defined for the most general type of object: anything that can be
broken down to atoms, i.e. atoms, molecules, collections, universes,
etc., i.e. in the class :class:`MMTK.Collections.GroupOfAtoms`.

The most elementary operations are inquiries about specific properties
of an object: number of atoms, total mass, center of mass, total momentum,
total charge, etc. There are also operations that compare two different
conformations of a system. Finally, there are special operations
for analyzing conformations of peptide chains and proteins.

Geometrical operations in periodic universes require special care.
Whenever a distance vector between two points in a systems is
evaluated, the minimum-image convention must be used in order to
obtain consistent results. MMTK provides routines for finding
these distance vectors as well as distances, angles, and dihedral
angles between any points. Because these operations depend on the
topology and geometry of the universe, they are implemented as
methods in class :class:`MMTK.Universe.Universe`
and its subclasses. Of course they are available for non-periodic
universes as well.

Universes also provide methods for obtaining :ref:`atom property
<atom_property>` objects that describe the state of the system
(configurations, velocities, masses), and for restoring the system
state from a :ref:`trajectory` file.

Energy evaluation
=================

Energy evaluation requires a force field, and therefore all the
methods in this section are defined only for universe objects, i.e. in
class :class:`MMTK.Universe.Universe`.  However, they all take an
optional arguments (anything that can be broken down into atoms) that
indicates for which subset of the universe the energy is to be
evaluated. In addition to the potential energy, energy gradients and
second derivatives (force constants) can be obtained, if the force
field implements them. There is also a method that returns a
dictionary containing the values for all the individual force field
terms, which is often useful for analysis.

Surfaces and volumes
====================

Surfaces and volumes can be analyzed for anything consisting of
atoms. Both quantities are defined by assigning a radius to each atom;
the surface of the resulting conglomerate of overlapping spheres is
taken to be the surface of the atom group. Atom radii for surface
determination are usually called "van der Waals radii", but there is
no unique method for determining them. MMTK uses the values from
:ref:`[Bondi1964] <Bondi1964>`. However, users can change these values for each
individual atom by assigning a new value to the attribute
``vdW_radius``.

The operations provided in :mod:`MMTK.MolecularSurface`
include basic surface and volume calculation, determination of
exposed atoms, and identification of contacts between two objects.

Miscellaneous operations
########################

Saving, loading, and copying objects
====================================

MMTK provides an easy way to store (almost) arbitrary objects in files
and retrieve them later. All objects of interest to users can be
stored, including chemical objects, collections, universes, normal
modes, configurations, etc. It is also possible to store standard
Python objects such as numbers, lists, dictionaries etc., as well as
practically any user-defined objects. Storage is based on the standard
Python module pickle.

Objects are saved with :func:`MMTK.save` and restored with
:func:`MMTK.load`.  If several objects are to be stored in a single
file, use tuples: ``save((object1, object2), filename)`` and
``object1, object2 = load(filename)`` to retrieve the objects.

Note that storing an object in a file implies storing all objects
referenced by it as well, such that the size of the file can become
larger than expected. For example, a configuration object contains
a reference to the universe for which it is defined. Therefore
storing a configuration object means storing the whole universe
as well. However, nothing is ever written twice to the same
file. If you store a list or a tuple containing a universe and
a configuration for it, the universe is written only once.

Frequently it is also useful to copy an object, such as a molecule or
a configuration. There are two functions (which are actually taken
from the Python standard library module copy) for this purpose, which
have a somewhat different behaviour for container-type objects (lists,
dictionaries, collections etc.). :func:`MMTK.copy` returns a copy of
the given object. For a container object, it returns a new container
object which contains the same objects as the original one. If the
intention is to get a container object which contains copies of the
original contents, then ``MMTK.deepcopy(object)`` should be used. For
objects that are not container-type objects, there is no difference
between the two functions.

Exporting to specific file formats and visualization
====================================================

MMTK can write objects in specific file formats that can be used by
other programs. Three file formats are supported: the PDB format,
widely used in computational chemistry, the DCD format for
trajectories, written by the programs CHARMM, X-Plor, and NAMDm, and
read by many visualization programs, and the VRML format, understood
by VRML browsers as a representation of a three-dimensional scene for
visualization. MMTK also provides a more general interface that can
generate graphics objects in any representation if a special module
for that representation exists. In addition to facilitating the
implementation of new graphics file formats, this approach also
permits the addition of custom graphics elements (lines, arrows,
spheres, etc.)  to molecular representations.

PDB, VRML, and DCD files
------------------------

Any chemical object, collection, or universe can be written to a PDB
or VRML file by calling the method ``writeToFile``, defined in class
:class:`MMTK.Collections.GroupOfAtoms`.  PDB files are read via the
class :class:`MMTK.PDB.PDBConfiguration`.  DCD files can be read by a
:class:`MMTK.DCD.DCDReader` object.  For writing DCD files, there is
the function :func:`MMTK.DCD.writeDCDPDB`, which also creates a
compatible PDB file without which the DCD file could not be
interpreted.

Special care must be taken to ensure a correct mapping of atom numbers
when reading from a DCD file. In MMTK, each atom object has a unique
identity and atom numbers, also used internally for efficiency, are
not strictly necessary and are not used anywhere in MMTK's application
programming interface. DCD file, however, simply list coordinates
sorted by atom number. For interpreting DCD files, another file must
be available which allows the identification of atoms from their
number and vice versa; this can for example be a PDB file.

When reading DCD files, MMTK assumes that the atom order in the DCD
file is identical to the internal atom numbering of the universe for
which the DCD file is read. This assumption is in general valid only
if the universe has been created from a PDB file that is compatible
with the DCD file, without any additions or removals.

Visualization and animation
---------------------------

The most common need for file export is visualization. All objects
that can be visualized (chemical systems and subsets thereof, normal
mode objects, trajectories) provide a method view
which creates temporary export files, starts a visualization program,
and deletes the temporary files. Depending on the object type there are
various optional parameters.

MMTK also allows visualization of normal modes and trajectories using
animation. Since not all visualization programs permit animation, and
since there is no standard way to ask for it, animation is implemented
only for the programs `XMol <http://www.msc.edu/msc/docs/xmol/>`_
and `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_. Animation is available for
normal modes, trajectories, and arbitrary sequences of configurations
(see function :func:`MMTK.Visualization.viewSequence`).

For more specialized needs, MMTK permits the creation of graphical
representations of most of its objects via general graphics modules that
have to be provided externally. Suitable modules are provided in the
package Scientific.Visualization and cover VRML (version 1), VRML2
(aka VRML97), and the molecular visualization program VMD. Modules for other
representations (e.g. rendering programs) can be written easily; it is
recommended to use the existing modules as an example. The generation
of graphics objects is handled by the method graphicsObjects,
defined in the class :class:`MMTK.Visualization.Viewable`,
which is a :term:`Mix-in class` that makes
graphics objects generation available for all objects that define
chemical systems or parts thereof, as well as for certain other objects
that are viewable.

The explicit generation of graphics objects permits the mixture of
different graphical representations for various parts of a system,
as well as the combination of MMTK-generated graphics objects with
arbitrary other graphics objects, such as lines, arrows, or spheres.
All graphics objects are finally combined into a scene object (also
defined in the various graphics modules) in order to be displayed.
See also the :ref:`Example-Visualization` examples.

Fields
======

For analyzing or visualizing atomic properties that change little over
short distances, it is often convenient to represent these properties as
functions of position instead of one value per atom. Functions of
position are also known as fields, and mathematical techniques for the
analysis of fields have proven useful in many branches of physics. Such
a field can be obtained by averaging over the values corresponding to
the atoms in a small region of space. MMTK provides classes for
scalar and vector field in module :mod:`MMTK.Field`.
See also the example :ref:`vector_field.py <Example-Miscellaneous>`.

Charge fitting
==============

A frequent problem in determining force field parameters is the
determination of partial charges for the atoms of a molecule by fitting
to the electrostatic potential around the molecule, which is obtained
from quantum chemistry programs. Although this is essentially a
straightforward linear least-squares problem, many procedures that are
in common use do not use state-of-the-art techniques and may yield
erroneous results. MMTK provides a charge fitting method that is
numerically stable and allows the imposition of constraints on the
charges. It is implemented in module :mod:`MMTK.ChargeFit`.
See also the example :ref:`charge_fit.py <Example-Miscellaneous>`.

.. _database:

Constructing the database
#########################

MMTK uses a database of chemical entities to define the properties of
atoms, molecules, and related objects. This database consists of plain
text files, more precisely short Python programs, whose names are the
names of the object types. This chapter explains how to construct and
manage these files. Note that the standard database already contains
many definitions, in particular for proteins and nucleic acids.
You do not need to read this chapter unless you want to add your
own molecule definitions.

MMTK's database does not have to reside in a single place. It can
consist of any number of subdatabases, each of which can be a
directory or a URL. Typically the database consists of at least two
parts: MMTK's standard definitions and a user's personal definitions.
When looking up an object type in the database, MMTK checks the value
of the environment variable MMTKDATABASE. The value of this variable
must be a list of subdatabase locations seperated by white space. If
the variable MMTKDATABASE is not defined, MMTK uses a default value
that contains the path ".mmtk/Database" in the user's home directory
followed by MMTK's standard database, which resides in the directory
Database within the MMTK package directory (on many Unix systems this
is ``/usr/local/lib/python2.x/site-packages/MMTK``).  MMTK checks the
subdatabases in the order in which they are mentioned in MMTKDATABASE.

Each subdatabase contains directories corresponding to the object
classes, i.e. Atoms (atom definitions), Groups (group definitions),
Molecules (molecule definitions), Complexes (complex definitions),
Proteins (protein definitions), and PDB (Protein Data Bank files).
These directories contain the definition files, whose names may
not contain any upper-case letters. These file names correspond
to the object types, e.g. the call ``MMTK.Molecule('Water')``
will cause MMTK to look for the file Molecules/water in the database
(note that the names are converted to lower case).

The remaining sections of this chapter explain how the individual
definition files are constructed. Keep in mind that each file is
actually a Python program, so of course standard Python syntax rules
apply.

Atom definitions
================

An atom definition in MMTK describes a chemical element, such as
"hydrogen". This should not be confused with the "atom types" used in
force field descriptions and in some modelling programs. As a
consequence, it is rarely necessary to add atom definitions to MMTK.

Atom definition files are short and of essentially identical format.
This is the definition for carbon:

::

    name = 'carbon'
    symbol = 'C'
    mass = [(12, 98.90), (13.003354826, 1.10)]
    color = 'black'
    vdW_radius = 0.17

The name should be meaningful to users, but is not used by MMTK
itself. The symbol, however, is used to identify chemical elements. It
must be exactly equal to the symbol defined by IUPAC, including
capitalization (e.g. 'Cl' for chlorine). The mass can be either a number
or a list of tuples, as shown above. Each tuple defines an isotope by
its mass and its percentage of occurrence; the percentages must add up
to 100. The color is used for VRML output and must equal one of the
color names defined in the module VRML. The van der Waals radius is used
for the calculation of molecular volumes and surfaces; the values are
taken from :ref:`[Bondi1964] <Bondi1964>`.

An application program can create an isolated atom with ``Atom('c')``
or, specifying an initial position, with ``Atom('c',
position=Vector(0.,1.,0.))``. The element name can use any combination
of upper and lower case letters, which are considered equivalent.

Group definitions
=================

Group definitions in MMTK exist to facilitate the definition of
molecules by avoiding the frequent repetition of common combinations.
MMTK doesn't give any physical meaning to groups. Groups can contain
atoms and other groups. Their definitions look exactly like molecule
definitions; the only difference between groups and molecules is the way
they are used.

This is the definition of a methyl group:

::

    name = 'methyl group'
    C  = Atom('C')
    H1 = Atom('H')
    H2 = Atom('H')
    H3 = Atom('H')
    bonds = [Bond(C, H1), Bond(C, H2), Bond(C, H3)]
    pdbmap = [('MTH', {'C': C, 'H1': H1, 'H2': H2, 'H3': H3})]
    amber_atom_type = {C: 'CT', H1: 'HC', H2: 'HC', H3: 'HC'}
    amber_charge = {C: 0., H1: 0.1, H2: 0.1, H3: 0.1}

The name should be meaningful to users, but is not used by MMTK
itself. The following lines create the atoms in the group and assign
them to variables. These variables become attributes of whatever
object uses this group; their names can be anything that is a legal
Python name. The list of bonds, however, must be assigned to the name
``bonds``. The bond list is used by force fields and for
visualization.

The name ``pdbmap`` is used for reading and writing PDB files. Its
value must be a list of tuples, where each tuple defines one PDB
residue. The first element of the tuple is the residue name, which is
used only for output. The second element is a dictionary that maps PDB
atom names to the actual atoms. The ``pdbmap`` entry of any object can be
overridden by an entry in a higher-level object. Therefore the entry for
a group is only used for atoms that do not occur in the entry for a
molecule that contains this group.

The remaining lines in the definition file contain information
specific to force fields, in this case the Amber force field. The
dictionary ``amber_atom_type`` defines the atom type for each atom;
the dictionary ``amber_charge`` defines the partial charges. As for
``pdbmap`` entries, these definitions can be overridden by
higher-level definitions.

Molecule definitions
====================

Molecules are typically used directly in application programs, but they
can also be used in the definition of complexes. Molecule definitions
can use atoms and groups.

This is the definition of a water molecule:

::

    name = 'water'
    structure = \
    "  O   \n" + \
    " / \  \n" + \
    "H   H \n"
    O  = Atom('O')
    H1 = Atom('H')
    H2 = Atom('H')
    bonds = [Bond(O, H1), Bond(O, H2)]
    pdbmap = [('HOH', {'O': O, 'H1': H1, 'H2': H2})]
    pdb_alternative = {'OH2': 'O'}
    amber_atom_type = {O: 'OW', H1: 'HW', H2: 'HW'}
    amber_charge = {O: -0.83400, H1: 0.41700, H2: 0.41700}
    configurations = {
    'default': ZMatrix([[H1], \
                       [O,  H1,  0.9572*Ang], \
                       [H2, O,   0.9572*Ang,  H1,  104.52*deg]])
    }

The name should be meaningful to users, but is not used by MMTK
itself. The structure is optional and not used by MMTK either. The
following lines create the atoms in the group and assign them to
variables. These variables become attributes of the molecule,
i.e. when a water molecule is created in an application program by ``w
= Molecule('water')``, then ``w.H1`` will refer to its first hydrogen
atom. The names of these variables can be any legal Python names. The
list of bonds, however, must be assigned to the name ``bonds``. The
bond list is used by force fields and for visualization.

The name ``pdbmap`` is used for reading and writing PDB files. Its
value must be a list of tuples, where each tuple defines one PDB
residue. The first element of the tuple is the residue name, which is
used only for output. The second element is a dictionary that maps PDB
atom names to the actual atoms. The ``pdbmap`` entry of any object can be
overridden by an entry in a higher-level object, i.e. in the case of a
molecule a complex containing it. The name ``pdb_alternative`` allows
to read PDB files that use non-standard names. When a
PDB atom name is not found in the ``pdbmap``, an attempt is made to
translate it to another name using ``pdb_alternative``.

The two following lines in the definition file contain information
specific to force fields, in this case the Amber force field. The
dictionary ``amber_atom_type`` defines the atom type for each atom; the
dictionary ``amber_charge`` defines the partial charges. As for pdbmap
entries, these definitions can be overridden by higher-level
definitions.

The name ``configurations`` can be defined to be a dictionary of
configurations for the molecule. During the construction of a molecule,
a configuration can be specified via an optional parameter, e.g.
``w = Molecule('water', configuration='default')``. The names of the
configurations can be arbitrary; only the name "default" has a special
meaning; it is applied by default if no other configuration is specified
when constructing the molecule. If there is no default configuration,
and no other configuration is explicitly specified, then the molecule is
created with undefined atomic positions.

There are three ways of describing configurations:

- By a Z-Matrix:

  ::

      ZMatrix([[H1], \
               [O,  H1,  0.9572*Ang], \
               [H2, O,   0.9572*Ang,  H1,  104.52*deg]])

- By Cartesian coordinates:

  ::

      Cartesian({O:  ( 0.004, -0.00518, 0.0),
      H1: (-0.092, -0.00518, 0.0),
      H2: ( 0.028,  0.0875,  0.0)})

- By a PDB file:

  ::

      PDBFile('water.pdb')

  The PDB file must be in the database subdirectory PDB, unless a full
  path name is specified for it.

Complex definitions
===================

Complexes are defined much like molecules, except that they are composed
of molecules and atoms; no groups are allowed, and neither are bonds.

Protein definitions
===================

Protein definitions can take many different forms, depending on the
source of input data and the type of information that is to be stored.
For proteins it is particularly useful that database definition files
are Python programs with all their flexibility.

The most common way of constructing a protein is from a PDB file. This
is an example for a protein definition:

::

    name = 'insulin'
    # Read the PDB file.
    conf = PDBConfiguration('insulin.pdb')
    # Construct the peptide chains.
    chains = conf.createPeptideChains()
    # Clean up
    del conf

The name should be meaningful to users, but is not used by MMTK
itself. The second command reads the sequences of all peptide chains
from a PDB file. Everything which is not a peptide chain is ignored.
The following line constructs a PeptideChain object (a special
molecule) for each chain from the PDB sequence. This involves
constructing positions for any missing hydrogen atoms.
Finally, the temporary data ("conf") is deleted, otherwise
it would remain in memory forever.

The net result of a protein definition file is the assignment of a list
of molecules (usually PeptideChain objects) to the variable "chains".
MMTK then constructs a protein object from it. To use the above example,
an application program would use the command
``p = Protein('insulin')``. The construction of the protein involves
one nontrivial (but automatic) step: the construction of disulfide
bridges for pairs of cystein residues whose sulfur atoms have a distance
of less then 2.5 Ångström.

Threads and parallelization
###########################

This chapter explains the use of threads by MMTK and MMTK's
parallelization support. This is an advanced topic, and not essential
for the majority MMTK applications. You need to read this chapter only
if you use multiprocessor computers, or if you want to implement
multi-threaded programs that use MMTK.

Threads are different execution paths through a program that are
executed in parallel, at least in principle; real parallel execution
is possible only on multiprocessor systems. MMTK makes use of threads
in two ways, which are conceptually unrelated: parallelization of
energy evaluation on shared-memory multiprocessor computers, and
support for multithreaded applications. Thread support is not
available on all machines; you can check if yous system supports
threads by starting a Python interpreter and typing import
threading. If this produces an error message, then your
system does not support threads, otherwise it is available in Python
and also in MMTK. If you do not have thread support in Python although
you know that your operating system supports threads, you might have
compiled your Python interpreter without thread support; in that case,
MMTK does not have thread support either.

Another approach to parallelization is message passing: several
processors work on a program and communicate via a fast network to
share results. A standard library, called MPI (Message Passing
Interface), has been developped for sharing data by message passing,
and implementations are available for all parallel computers currently
on the market. MMTK contains elementary support for parallelization by
message passing: only the energy evaluation has been paralellized,
using a data-replication strategy, which is simple but not the most
efficient for large systems. MPI support is disabled by default.
Enabling it involves modifying the file Src/Setup.template prior
to compilation of MMTK. Furthermore, an MPI-enabled installation of
ScientificPython is required, and the mpipython executable must
be used instead of the standard Python interpreter.

Threads and message passing can be used together to use a cluster of
shared-memory machines most efficiently. However, this requires that
the thread and MPI implementations being used work together; sometimes
there are conflicts, for example due to the use of the same signal in
both libraries. Refer to your system documentation for details.

The use of threads for parallelization on shared-memory systems is
very simple: Just set the environment variable MMTK_ENERGY_THREADS to
the desired value.  If this variable is not defined, the default value
is 1, i.e. energy evaluations are performed serially. For choosing an
appropriate value for this environment variable, the following points
should be considered:

- The number of energy evaluation threads should not be larger than the
  number of processors that are fully dedicated to the MMTK application.
  A larger number of threads does not lead to wrong results,
  but it can increase the total execution time.

- MMTK assumes that all processors are equally fast. If you use a
  heteregenous multiprocessor machine, in which the processors have
  different speeds, you might find that the total execution time is
  larger than without threads.

- The use of threads incurs some computational overhead. For very small
  systems, it is usually faster not to use threads.

- Not all energy terms necessarily support threads. Of the force field
  terms that part of MMTK, only the multipole algorithms for
  electrostatic interactions does not support threads, but additional
  force fields defined outside MMTK might also be affected. MMTK
  automatically evaluates such energy terms with a single thread, such
  that there is no risk of getting wrong results. However, you might not
  get the performance you expect.

- If second derivatives of the potential energy are requested, energy
  evaluation is handled by a single thread. An efficient implementation
  of multi-threaded energy evaluation would require a separate copy of
  the second-derivative matrix per thread. This approach needs too much
  memory for big systems to be feasible. Since second derivatives are
  almost exclusively used for normal mode calculations, which need only
  a single energy evaluation, multi-thread support is not particularly
  important anyway.

Parallelization via message passing is somewhat more complicated.
In the current MMTK parallelization model, all processors execute
the same program and replicate all tasks, with the important exception
of energy evaluation. Energy terms are divided evenly between the
processors, and at the end the energy and gradient values are shared
by all machines. This is the only step involving network communication.
Like thread-based parallelization, message-passing parallelization
does not support the evaluation of second derivatives.

A special problem with message-passing systems is input and output.
The MMTK application must ensure that output files are written by
only one processor, and that all processors correctly access input
files, especially in the case of each processor having its own
disk space. See the example :ref:`md.py <Example-MPI>`
for illustration.

Multithreaded applications are applications that use multiple threads
in order to simplify the implementation of certain algorithms, i.e.
not necessarily with the goal of profiting from multiple processors.
If you plan to write a multithreaded application that uses MMTK,
you should first make sure you understand threading support in
Python. In particular, you should keep in mind that the global
interpreter lock prevents the effective use of multiple processors
by Python code; only one thread at a time can execute interpreted
Python code. C code called from Python can permit other threads
to execute simultaneously; MMTK does this for energy evaluation,
molecular dynamics integration, energy minimization, and normal
mode calculation.

A general problem in multithreaded applications is access to resources
that are shared among the threads. In MMTK applications, the most
important shared resource is the description of the chemical systems,
i.e. universe objects and their contents. Chaos would result if two
threads tried to modify the state of a universe simultaneously, or
even if one thread uses information that is simultaneously being
modified by another thread. Synchronization is therefore a critical
part of multithreaded application. MMTK provides two synchronization
aids, both of which described in the documentation of the class 
:class:`MMTK.Universe.Universe`: the configuration change
lock (methods :func:`~MMTK.Universe.Universe.acquireConfigurationChangeLock`
and :func:`~MMTK.Universe.Universe.releaseConfigurationChangeLock`), 
and the universe state lock (methods 
:func:`~MMTK.Universe.Universe.acquireReadStateChangeLock`, 
:func:`~MMTK.Universe.Universe.releaseReadStateChangeLock`,
:func:`~MMTK.Universe.Universe.acquireWriteStateChangeLock`, and
:func:`~MMTK.Universe.Universe.releaseWriteStateChangeLock`). 
Only a few common universe operations manipulate the universe 
state lock in order to avoid conflicts with other threads; 
these methods are marked as thread-safe in the description. 
All other operations should only be
used inside a code section that is protected by the appropriate
manipulation of the state lock. The configuration change lock is less
critical; it is used only by the molecular dynamics and energy
minimization algorithms in MMTK.

Bibliography
############

.. _Bondi1964:

[Bondi1964]

  | A. Bondi
  | `van der Waals Volumes and Radii <http://dx.doi.org/10.1021/j100785a001>`_
  | 1964

.. _Eisenhaber1993:

[Eisenhaber1993]

  | F. Eisenhaber, P. Argos
  | `Improved Strategy in Analytic Surface Calculation for Molecular
    Systems: Handling of Singularities and Computational Efficiency
    <http://dx.doi.org/10.1002/jcc.540141103>`_
  | 1993

.. _Eisenhaber1995:

[Eisenhaber1995]

 | F. Eisenhaber, P. Lijnzaad, P. Argos, M. Scharf
 | `The Double Cubic Lattice Method: Efficient Approaches to Numerical
    Integration of Surface Area and Volume and to Dot Surface
    Contouring of Molecular Assemblies
    <http://dx.doi.org/10.1002/jcc.540160303(YoYo)>`_
 | 1995

.. _Hinsen1995:

[Hinsen1995]

 | Konrad Hinsen, Gerald R. Kneller
 | `Influence of constraints on the dynamics of polypeptide chains
    <http://dx.doi.org/10.1103/PhysRevE.52.6868>`_
 | 1995

.. _Hinsen1997:

[Hinsen1997]

 | Konrad Hinsen, Benoit Roux
 | `An accurate potential for simulating proton transfer in acetylacetone
    <http://dx.doi.org/10.1002/(SICI)1096-987X(199702)18:3%3C368::AID-JCC7%3E3.0.CO;2-S>`_
 | 1997

.. _Hinsen1998:

[Hinsen1998]

 | Konrad Hinsen
 | `Analysis of domain motions by approximate normal mode calculations
    <http://dx.doi.org/10.1002/(SICI)1097-0134(19981115)33:3%3C417::AID-PROT10%3E3.0.CO;2-8>`_
 | 1998

.. _Hinsen1999:

[Hinsen1999]

 | Konrad Hinsen, Aline Thomas, Martin J. Field
 | `Analysis of domain motions in large proteins
    <http://dx.doi.org/10.1002/(SICI)1097-0134(19990215)34:3%3C369::AID-PROT9%3E3.0.CO;2-F>`_
 | 1999

.. _Hinsen1999a:

[Hinsen1999a]

 | Konrad Hinsen, Gerald R. Kneller
 | `Projection methods for the analysis of complex motions in macromolecules
    <http://dx.doi.org/10.1080/08927020008025373>`_
 | 1999

.. _Hinsen1999b:

[Hinsen1999b]

 | Konrad Hinsen, Gerald R. Kneller
 | `A simplified force field for describing vibrational protein
    dynamics over the whole frequency range
    <http://dx.doi.org/10.1063/1.480441>`_
 | 1999

.. _Hinsen2000:

[Hinsen2000]

 | Konrad Hinsen, Andrei J. Petrescu, Serge Dellerue,
   Marie-Claire Bellissent-Funel, Gerald R. Kneller.
 |  `Harmonicity in slow protein dynamics
     <http://dx.doi.org/10.1016/S0301-0104(00)00222-6>`_
 | 2000

.. _Kneller1990:

[Kneller1990]

 | Gerald R. Kneller
 | `Superposition of molecular structures using quaternions
    <http://dx.doi.org/10.1080/08927029108022453>`_
 | 1990

.. _Kneller1996:

[Kneller1996]

 | Gerald R. Kneller, Thomas Mülders
 | `Nosé-Andersen dynamics of partially rigid molecules:
    Coupling of all degrees of freedom to heat and pressure baths
    <http://dx.doi.org/10.1103/PhysRevE.54.6825>`_
 | 1996

.. _Swope1982:

[Swope1982]

 | W.C. Swope, H.C. Andersen, P.H. Berens, K.R. Wilson
 | `A computer simulation method for the calculation of equilibrium 
    constants for the formation of physical clusters of molecules:
    application to small water clusters
    <http://dx.doi.org/10.1063/1.442716>`_
 | 1982

.. _vanGunsteren1982:

[vanGunsteren1982]

 | Wilfred F. van Gunsteren, Martin Karplus
 | `Effect of Constraints on the Dynamics of Macromolecules
    <http://dx.doi.org/10.1021/ma00234a015>`_
 | 1982

.. _Viduna2000:

[Viduna2000]

 | David Viduna, Konrad Hinsen, Gerald R. Kneller
 |  `The influence of molecular flexibility on DNA radiosensitivity:
     A simulation study <http://dx.doi.org/10.1103/PhysRevE.62.3986>`_
 | 2000

.. _Wolf1999:

[Wolf1999]

 | D. Wolf, P. Keblinski, S.R. Philpot, J. Eggebrecht
 | `Exact method for the simulation of Coulombic systems by spherically
    truncated, pairwise r\ :sup:`-1` summation
    <http://dx.doi.org/10.1063/1.478738>`_
 | 1999
