This repository provides the SedFoam solver and the new extension, overSedDymFoam, which handles dynamic objects and incorporate the 6DoF in sedFoam.

Features
--------
A three-dimensional two-phase flow solver, SedFoam, has been developed for sediment transport applications. The solver is extended from twoPhaseEulerFoam available in the 2.1.0 release of the open-source CFD (computational fluid dynamics) toolbox OpenFOAM. Currently the latest version is SedFoam-2306. In this approach the sediment phase is modeled as a continuum, and constitutive laws have to be prescribed for the sediment stresses. Within this solver, we have integrated two distinct intergranular stress models: the Kinetic Theory of Granular Flows and the Dense Granular Flow Rheology μ(I). When it comes to modeling fluid stresses, it has the versatility to simulate laminar or turbulent flow regimes. Furthermore, there are three turbulence models at your disposal for sediment transport simulations, including a straightforward mixing length model (restricted to one-dimensional configurations), a k-ε model, the classical k-ω model, and the k-ω model by Wilcox (2006). In addition to these features, the extension known as overSedDymFoam equips SedFoam with the capability to account for the dynamic movement of objects due to both hydrodynamic and granular stresses. To delve deeper into this extension, you can refer to the article titled "The Role of Moving Objects in Sediment Transport Processes" in the OpenFOAM journal. The capabilities of overSedDymFoam are integrated into the following directories, alongside with their respective files, which are now part of the SedFoam ecosystem:

- `forcesSedFoam`: this folder allows the computation of hydrodynamic and granular stresses exerted on solid objects.
- `sixDoFRigidBodyMotionSedFoam`: this folder contains the 6dof solver that allows the motion of a rigid body based on the computed forces acting on it.
- `solver/overSedDymFoam`: this is the sedFoam solver modified to deal with dynamic meshes.

Installation
------------
```bash
cd sedfoam
./Allwclean
./Allwmake
```

Data reproducibility
-----
Numerical cases presented in this article are located in `OF_article_cases`. Numerical simulations are performed using OF v2306 release and can be easily executed running `./Allrun` for each case. All cases are set to run in parallel and the number or processors can be reduced if needed. The following cases are ready to reproduce the numerical results of the manuscript:
- `FallingSphereOverset` and `FallingSphereMorphing` correspond to the falling sphere in a viscous fluid using an overset and a morphing mesh respectively.
- `FallingSphereSuspensionOverset` and `FallingSphereSuspensionMorphing` correspond to the falling sphere in a viscous fluid with a dilute particle suspension using an overset and a morphing mesh respectively.
- `RestingSphereMorphing` correspond to the falling sphere on a horizontal compacted granular bed using a morphing mesh.
- `2DCylinderUniformGranularFlow` corresponds to the 2D cylinder immersed in a uniform granular flow experiencing lift using a morphing mesh.
Postprocessing plots can be generated using python. To do so, you just have to run the python scripts located in the folder `OF_article_cases/Py`.

Additional tutorials can be found in `tutorialsDyM`. This folder contains several examples with moving objects interacting with sediment transport processes but not all of them are part of the manuscript "The Role of Moving Objects in Sediment Transport Processes" presented in the OpenFOAM journal.

It is worth mentioning that the python package fluidfoam ( https://fluidfoam.readthedocs.io/ ) is needed for postprocessing of the tutorials. One can easily install fluidfoam by executing the following command:

```bash
pip install fluidfoam --user
```
