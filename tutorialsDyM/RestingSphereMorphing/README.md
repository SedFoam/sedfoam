 Resting sphere - Morphing approach
 ============

 This tutorial tests the ability of overSedDymFoam to capture the particle pressure and frictional stresses on a rigid body. We consider a rigid sphere with a diameter of D=0.015m  falling into a dense granular bed. The computational domain is assumed to be spherical. A horizontal sediment bed is included in the numerical simulation. The bed interface is located at a distance of 0.84D from the lowest point of the sphere.


Usage
-----

To speed up the numerical simulation, we need to run the initial_1D_profile/1DSedim case to get the sedimentation of the granular material. Once we are at the tutorialsDyM/RestingSphereMorphing/initial_1D_profile/1DSedim directory, we should execute the following command:
```bash
./Allrun
```
Once the initial_1D_profile/1DSedim numerical simulation is over, we should come back to tutorialsDyM/RestingSphereMorphing folder and execute the following command to launch the 3D numerical case:

```bash
./Allrun
```

Postprocessing
---------
To generate the postprocessing files you need to execute makeFiles in the case first. 

```bash
./makeFiles
```

Then, we can run the python script plotRestingSphere.py located in the folder tutorialsDyM/Py.

```bash
python plotRestingSphere.py
```
