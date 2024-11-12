sedfoam
=======

[![Release](https://img.shields.io/badge/release-2406-blue.svg)](http://github.com/SedFoam/sedfoam)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/b9ad60ec6171496290c336697426cd48)](https://www.codacy.com/gh/SedFoam/sedfoam/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=SedFoam/sedfoam&amp;utm_campaign=Badge_Grade)
[![OpenFOAM v24xx](https://img.shields.io/badge/OpenFOAM-v24xx-brightgreen.svg)](https://openfoam.com/)
[![OpenFOAM v23xx](https://img.shields.io/badge/OpenFOAM-v23xx-brightgreen.svg)](https://openfoam.com/)
[![OpenFOAM v22xx](https://img.shields.io/badge/OpenFOAM-v22xx-brightgreen.svg)](https://openfoam.com/)
[![OpenFOAM v21xx](https://img.shields.io/badge/OpenFOAM-v21xx-brightgreen.svg)](https://openfoam.com/)
[![OpenFOAM v20xx](https://img.shields.io/badge/OpenFOAM-v20xx-brightgreen.svg)](https://openfoam.com/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.836642.svg)](https://doi.org/10.5281/zenodo.836642)

This repository provides the SedFoam solver.

[![](https://i.ibb.co/WgS6PYB/Capture-d-e-cran-2018-12-14-a-11-27-01.png)](https://www.youtube.com/watch?v=cVf7qm_ZDK0)

Status
------

The sedfoam solver is in development and is not yet fully functional.

Pull requests are encouraged!

Features
--------
A three-dimensional two-phase flow solver, SedFoam, has been developed for sediment transport applications. The solver is extended from twoPhaseEulerFoam available in the 2.1.0 release of the open-source CFD (computational fluid dynamics) toolbox OpenFOAM. In this approach the sediment phase is modeled as a continuum, and constitutive laws have to be prescribed for the sediment stresses. In the proposed solver, two different intergranular stress models are implemented: the kinetic theory of granular flows and the dense granular flow rheology μ(I). For the fluid stress, laminar or turbulent flow regimes can be simulated and three different turbulence models are available for sediment transport: a simple mixing length model (one-dimensional configuration only), a k − ε, a classical k − ω and a k − ω from Wilcox (2006).

Installation
------------

```bash
cd $WM_PROJECT_USER_DIR
git clone --recurse-submodules https://github.com/sedfoam/sedfoam sedfoam
cd sedfoam
./Allwclean
./Allwmake
```

Usage
-----

There are tutorials located in `tutorials` and details can be found [here.](http://sedfoam.github.io/sedfoam)

How to cite
-----------

The sedfoam solver should be referenced using the GMD paper :

Julien Chauchat, Zhen Cheng, Tim Nagel, Cyrille Bonamy, and Tian-Jian Hsu (2017) SedFoam-2.0: a 3-D two-phase flow numerical model for sediment transport (Geosci. Model Dev., 10, 4367-4392) [![DOI](https://img.shields.io/badge/DOI-10.5195%2Fgmd_10_4367_2017-blue.svg)](https://doi.org/10.5194/gmd-10-4367-2017)

and the reference in which a specific configuration and/or feature has been first presented (please see list of references in the documentation).

If you are using a sedFoam version without modification, you should also refer to the Zenodo DOI of the version, current version is [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.836642.svg)](https://doi.org/10.5281/zenodo.836642)


Publications
------------

Julien Chauchat, Zhen Cheng, Tim Nagel, Cyrille Bonamy, and Tian-Jian Hsu (2017) SedFoam-2.0: a 3-D two-phase flow numerical model for sediment transport (Geosci. Model Dev., 10, 4367-4392) [![DOI](https://img.shields.io/badge/DOI-10.5195%2Fgmd_10_4367_2017-blue.svg)](https://doi.org/10.5194/gmd-10-4367-2017)

Antoine Mathieu, Julien Chauchat, Cyrille Bonamy, Tim Nagel (2019) Two-phase flow simulation of tunnel and lee-wake erosion of scour below a submarine pipeline (Water, 11, 8) [![DOI](https://img.shields.io/badge/DOI-10.3390%2Fw11081727-blue.svg)](https://www.mdpi.com/2073-4441/11/8/1727)

Tim Nagel, Julien Chauchat, Cyrille Bonamy, Xiaofeng Liu, Zhen Cheng, and Tian-Jian Hsu (2020) Three-dimensional scour simulations with a two-phase flow model (Advances in Water Resources) [![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.advwatres.2020.103544-blue.svg)](https://doi.org/10.1016/j.advwatres.2020.103544) 

Developers
----------

*   Cyrille Bonamy
*   Julien Chauchat
*   Antoine Mathieu
*   Eduard Puig Montellà
*   Rémi Chassagne
*   Zhen Cheng
*   Tim Nagel
*   Benjamin Tsai
*   Tian-Jian Hsu

Contributors
------------

*   Ali Salimi-Tarazouj
*   Hugo Divel
*   Pablo Higuera

Acknowledgements
----------------

[OpenFOAM](https://www.openfoam.com) is the free, open source CFD
software developed primarily by OpenCFD Ltd since 2004.
The OpenFOAM trademark is owned by OpenCFD Ltd.

![scour image](https://i.ibb.co/pWjZqd4/scour3-D-cylinder.jpg)
