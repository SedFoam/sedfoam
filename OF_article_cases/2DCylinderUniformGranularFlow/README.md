Uniform granular flow around a cylinder - Morphing approach
 ============

In this tutorial we analyze  the trajectory of an intruder immersed in a uniform granular flow. To reduce the computational cost, we use a cylinder (a 2D object in OpenFOAM). This tutorial is inspired by the work of:

Guillard, François, Yoël Forterre, and Olivier Pouliquen. "Lift forces in granular media." Physics of Fluids 26.4 (2014),

where lift and drag forces on a cylinder were measured under a granular flow. We adopt a similar setup: a granular packing is considered, then, the bottom plate is set to a prescribed constant velocity of v=0.01m/s that induces a uniform granular flow above the bottom plate. Placing a cylinder in the granular medium distorts the flow patterns and generates a lift force on the cylinder. The initial velocity field is set to v=0.01m/s as well.



Usage
-----

To speed up the numerical simulation, we need to run the Case1D to get the sedimentation of the granular material. Once we are at the tutorialsDyM/2DCylinderUniformGranularFlow/Case1D directory, we should execute the following command:
```bash
./Allrun
```
Once the Case1D numerical simulation is over, we should come back to tutorialsDyM/2DCylinderUniformGranularFlow folder and execute the following command to launch the 2D numerical case:

```bash
./Allrun
```

Postprocessing
---------

We can run the python script plot2DCylinderPosition.py located in the folder tutorialsDyM/Py.

```bash
python plot2DCylinderPosition.py
```
