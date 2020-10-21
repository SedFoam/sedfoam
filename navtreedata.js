var NAVTREE =
[
  [ "SedFOAM", "index.html", [
    [ "Installation", "install.html", [
      [ "Compiling from source with OpenFOAM-5.0", "install.html#install_src_unix", null ],
      [ "Compiling from source with FOAM-extend release 4.0", "install.html#install_src_unix_foamext", null ]
    ] ],
    [ "Getting started", "starting.html", [
      [ "Running a case", "starting.html#running", null ]
    ] ],
    [ "Tutorials", "tutorials.html", [
      [ "1D test cases", "tutorials.html#One-dimensional", [
        [ "1DSedim: Pure sedimentation", "tutorials.html#Sedim_testcase", [
          [ "Pre-processing", "tutorials.html#preproc", [
            [ "Mesh generation", "tutorials.html#meshgen", null ],
            [ "Boundary and initial conditions", "tutorials.html#BC_init", null ],
            [ "Physical properties", "tutorials.html#PhysProp", null ],
            [ "Control", "tutorials.html#control_sedFoam", null ],
            [ "Discretisation and linear-solver settings", "tutorials.html#Disc_solv", null ]
          ] ],
          [ "Running sedFoam_rbgh", "tutorials.html#run_sedFoam", null ],
          [ "Post-processing", "tutorials.html#post", [
            [ "Using paraview", "tutorials.html#postparaview", null ],
            [ "Using fluidfoam", "tutorials.html#postfluidfoam", null ]
          ] ]
        ] ],
        [ "1DBedload: Laminar bedLoad", "tutorials.html#BedLoadCLB", [
          [ "Pre-processing", "tutorials.html#preprocBL", [
            [ "Mesh generation", "tutorials.html#meshgenBL", null ],
            [ "Boundary and initial conditions", "tutorials.html#BC_initiBL", null ],
            [ "Physical properties", "tutorials.html#PhysPropBL", null ],
            [ "Discretisation and linear-solver settings", "tutorials.html#Disc_solvBL", null ],
            [ "Control", "tutorials.html#control_sedFoamBL", null ]
          ] ],
          [ "Post-processing using python", "tutorials.html#postBL", null ],
          [ "Extension to the Mu(Iv) rheology", "tutorials.html#ExtensionBL", null ]
        ] ],
        [ "1DAvalancheMuI: Dry granular avalanche", "tutorials.html#dryavalanche", [
          [ "Pre-processing", "tutorials.html#preprocAva", [
            [ "Mesh generation", "tutorials.html#meshgenAva", null ],
            [ "Boundary and initial conditions", "tutorials.html#BC_initiAva", null ],
            [ "Physical properties", "tutorials.html#PhysPropAva", null ]
          ] ],
          [ "Post-processing using python", "tutorials.html#postAvalanche", null ]
        ] ],
        [ "1DSimpleShear: Simple shear with kinetic theory", "tutorials.html#SimpleShear", [
          [ "Post-processing using python", "tutorials.html#postSimpleShear", null ]
        ] ],
        [ "1DBoundaryLayer: Single-phase turbulent boundary layer flow", "tutorials.html#BoundaryLayer", null ],
        [ "1DSheetFlow: Turbulent sheet-flows", "tutorials.html#SheetFlow", [
          [ "Pre-processing", "tutorials.html#preprocSF", [
            [ "Mesh generation", "tutorials.html#meshgenSF", null ],
            [ "Boundary and initial conditions", "tutorials.html#BC_initiSF", null ],
            [ "Physical properties", "tutorials.html#PhysPropSF", null ],
            [ "Discretisation and linear-solver settings", "tutorials.html#Disc_solvSF1", null ],
            [ "Control", "tutorials.html#control_sedFoamSF", null ]
          ] ],
          [ "Post-processing using python", "tutorials.html#postSF", null ]
        ] ]
      ] ],
      [ "2D/3D test cases", "tutorials.html#Multi-dimensional", [
        [ "2DScour: Scour downstream of an apron", "tutorials.html#ScourApron", [
          [ "Pre-processing", "tutorials.html#preprocSA", null ],
          [ "Mesh generation", "tutorials.html#meshgenSA", null ],
          [ "Boundary and initial conditions", "tutorials.html#BC_initiSA", null ],
          [ "Physical properties", "tutorials.html#PhysPropSA", null ],
          [ "Discretisation and linear-solver settings", "tutorials.html#Disc_solvSF2", null ],
          [ "Control", "tutorials.html#control_sedFoamSF2", null ],
          [ "Computation launching", "tutorials.html#lauch_SA", null ],
          [ "Post-processing using python", "tutorials.html#postSF2", null ]
        ] ],
        [ "2DPipelineScour: Scour around an horizontal cylinder", "tutorials.html#Scour2DCylinder", null ],
        [ "3DScourCylinder: Scour around a vertical cylinder", "tutorials.html#Scour3DCylinder", null ]
      ] ]
    ] ],
    [ "Input description", "inputs.html", [
      [ "Physical processes", "inputs.html#inputParams", [
        [ "File transportProperties", "inputs.html#transportProperties", null ],
        [ "File interfacialProperties", "inputs.html#interfacialProperties", null ],
        [ "File ppProperties", "inputs.html#ppProperties", null ],
        [ "File granularRheologyProperties", "inputs.html#granularRheologyProperties", null ],
        [ "File forceProperties", "inputs.html#forceProperties", null ],
        [ "File turbulenceProperties", "inputs.html#turbulenceProperties", null ],
        [ "File twophaseRASProperties", "inputs.html#twophaseRASProperties", null ]
      ] ]
    ] ],
    [ "Download", "^http://github.com/sedfoam/sedfoam", null ]
  ] ]
];

var NAVTREEINDEX =
[
"index.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';