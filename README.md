sedfoam
=======

This repository provides the SedFoam solver.

Features
--------
A three-dimensional two-phase flow solver, SedFoam, has been developed for sediment transport applications. The solver is extended from twoPhaseEulerFoam available in the 2.1.0 release of the open-source CFD (computational fluid dynamics) toolbox OpenFOAM. In this approach the sediment phase is modeled as a continuum, and constitutive laws have to be prescribed for the sediment stresses. In the proposed solver, two different intergranular stress models are implemented: the kinetic theory of granular flows and the dense granular flow rheology μ(I). For the fluid stress, laminar or turbulent flow regimes can be simulated and three different turbulence models are available for sediment transport: a simple mixing length model (one-dimensional configuration only), a k − ε, a classical k − ω and a k − ω from Wilcox (2006). The last extension, overSedDymFoam, includes the ability of objects to move due to both hydrodynamic and granular stresses.

Installation
------------

```bash
cd sedfoam
./Allwclean
./Allwmake
```

Usage
-----

There are tutorials located in `tutorials`. Whereas the tutorials regarding dynamic objects are located in  `tutorialsDyM`.

