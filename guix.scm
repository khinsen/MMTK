(use-modules (guix)
             (guix build-system python)
             ((guix licenses) #:prefix license:)
             (guix git-download)
             (gnu packages python)
             (gnu packages python-xyz)
             (gnu packages maths))

(package
  (name "python2-mmtk")
  (version "2.7.12")
  (source (local-file (dirname (current-filename))
                      #:recursive? #t))
  (build-system python-build-system)
  (native-inputs
   `(("netcdf" ,netcdf)))
  (propagated-inputs
   `(("python-scientific" ,python2-scientific)
     ("python-tkinter" ,python-2 "tk")))
  (arguments
   `(#:python ,python-2
     #:tests? #f
     #:phases
     (modify-phases %standard-phases
       (add-before 'build 'includes-from-scientific
         (lambda* (#:key inputs #:allow-other-keys)
           (mkdir-p "Include/Scientific")
           (copy-recursively
            (string-append
             (assoc-ref inputs "python-scientific")
             "/include/python2.7/Scientific")
            "Include/Scientific"))))))
  (home-page "http://dirac.cnrs-orleans.fr/MMTK")
  (synopsis "Python library for molecular simulation")
  (description "MMTK is a library for molecular simulations with an emphasis
on biomolecules.  It provides widely used methods such as Molecular Dynamics
and normal mode analysis, but also basic routines for implementing new methods
for simulation and analysis.  The library is currently not actively maintained
and works only with Python 2 and NumPy < 1.9.")
  (license license:cecill-c))
