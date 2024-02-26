 Falling sphere in a sediment suspension - Morphing approach
 ============

This tutorial investigates the free fall of a sphere with a diameter of 15mm settling in silicon oil with  with the presence of sediments in suspension using the morphing mesh approach. For this numerical setup, the sediment content is set to 5%, and the particles are neutrally buoyant (i.e. same density as the fluid) with a mean diameter of d=0.29mm.


Usage
-----

The numerical case can be executed running the following command:
```bash
./Allrun
```

Postprocessing
---------

You can run the python script plot3DFallingSphere_suspension.py located in the folder tutorialsDyM/Py. Make sure that the chosen mesh in the script is set to DynamycMesh="Morphing".

```bash
python plot3DFallingSphere_suspension.py
```
