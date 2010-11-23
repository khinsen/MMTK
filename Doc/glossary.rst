Glossary
========

.. glossary::

   Abstract base class
      A :term:`Base Class` that is not
      directly usable by itself, but which defines the common properties
      of several subclasses. Example: the class 
      :class:`MMTK.ChemicalObjects.ChemicalObject` is
      an abstract base class which defines the common properties of its
      subclasses :class:`MMTK.ChemicalObjects.Atom`, 
      :class:`MMTK.ChemicalObjects.Group`, 
      :class:`MMTK.ChemicalObjects.Molecule`, 
      :class:`MMTK.ChemicaObjects.Complex`, and 
      :class:`MMTK.ChemicalObjects.AtomCluster`. A :term:`Mix-in class` 
      is a special kind of abstract base class.

   Base class
      A class from which another class inherits. In most
      cases, the inheriting class is a specialization of the base class.
      For example, the class :class:`MMTK.ChemicalObjects.Molecule` is a
      base class of :class:`MMTK.Proteins.PeptideChain`,
      because peptide chains are special molecules. Another common
      application is the :term:`Abstract base class`.

   Mix-in class
      A class that is used as a :term:`Base class`
      in other classes with the sole intention of providing methods
      that are common to these classes. Mix-in classes cannot be used
      to create instances. They are a special kind of
      :term:`Abstract base class`.
      Example: class :class:`MMTK.Collections.GroupOfAtoms`.

   Subclass
      A class that has another class as its :term:`Base class`.
      The subclass is usually a specialization of the base class, and can
      use all of the methods defined in the base class.
      Example: class :class:`MMTK.Proteins.Residue` is
      a subclass of :class:`MMTK.ChemicalObjects.Group`.


