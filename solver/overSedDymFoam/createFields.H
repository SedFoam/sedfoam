/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

IOdictionary twophaseRASProperties
(
    IOobject
    (
        "twophaseRASProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

IOdictionary ppProperties
(
    IOobject
    (
        "ppProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Read ppProperties


dimensionedScalar alphaMax
(
    ppProperties.lookupOrDefault
    (
        "alphaMax",
        dimensionedScalar("alphaMax",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.635)
    )
);

dimensionedScalar alphaMinFriction
(
    ppProperties.lookupOrDefault
    (
        "alphaMinFriction",
        dimensionedScalar("alphaMinFriction",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.57)
    )
);

dimensionedScalar Fr_
(
    ppProperties.lookupOrDefault
    (
        "Fr",
        dimensionedScalar("Fr",
                          dimensionSet(1, -1, -2, 0, 0, 0, 0),
                          5e-2)
    )
);
dimensionedScalar eta0_
(
    ppProperties.lookupOrDefault
    (
       "eta0",
       dimensionedScalar("eta0", dimensionSet(0, 0, 0, 0, 0, 0, 0), 3 )
    )
);
dimensionedScalar eta1_
(
    ppProperties.lookupOrDefault
    (
       "eta1",
       dimensionedScalar("eta1", dimensionSet(0, 0, 0, 0, 0, 0, 0), 5 )
    )
);
Switch packingLimiter(ppProperties.lookup("packingLimiter"));

// Read transportProperties

dimensionedScalar alphaDiffusion
(
    transportProperties.lookupOrDefault
    (
        "alphaDiffusion",
        dimensionedScalar("alphaDiffusion",
                          dimensionSet(0, 2, -1, 0, 0, 0, 0),
                          0e0)
    )
);

dimensionedScalar alphaSmall
(
    transportProperties.lookupOrDefault
    (
        "alphaSmall",
        dimensionedScalar("alphaSmall",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          1e-6)
    )
);

dimensionedScalar nuMax
(
    transportProperties.lookupOrDefault
    (
        "nuMax",
        dimensionedScalar("nuMax",
                          dimensionSet(0, 1, -2, 0, 0, 0, 0),
                          1e0)
    )
);

dimensionedScalar Dsmall
(
    transportProperties.lookupOrDefault
    (
        "Dsmall",
        dimensionedScalar("Dsmall",
                          dimensionSet(0, 0, -1, 0, 0, 0, 0),
                          1e-6)
    )
);

dimensionedScalar nuMaxExp
(
    transportProperties.lookupOrDefault
    (
        "nuMaxExp",
        dimensionedScalar("nuMaxExp",
                          dimensionSet(0, 1, -2, 0, 0, 0, 0),
                          1e0)
    )
);
// Initialize phase a and b variables

autoPtr<phaseModel> phasea = phaseModel::New
(
    mesh,
    transportProperties,
    "a"
);

autoPtr<phaseModel> phaseb = phaseModel::New
(
    mesh,
    transportProperties,
    "b"
);

volScalarField& alpha = phasea->alpha();
volVectorField& Ub = phaseb->U();
volScalarField& beta = phaseb->alpha();
Info<< "Warning : alpha_b read but not used!\n" << endl;
Info<< "alpha_b is set to 1 - alpha_a\n" << endl;
beta = 1 - alpha;
surfaceScalarField& phib = phaseb->phi();
const dimensionedScalar& rhob = phaseb->rho();
const dimensionedScalar& nub = phaseb->nu();

volVectorField& Ua = phasea->U();
surfaceScalarField& phia = phasea->phi();
const dimensionedScalar& rhoa = phasea->rho();
const dimensionedScalar& nua = phasea->nu();

// Reading ForceProperties

Info<< "\nReading forceProperties for gradP\n" <<endl;

IOdictionary forceProperties
(
    IOobject
    (
        "forceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar tilt
(
    forceProperties.lookupOrDefault
    (
        "tilt",
        dimensionedScalar("tilt", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
    )
);
dimensionedVector gradPAMP1
(
    dimensionedVector::getOrDefault
    (
        "gradPAMP1",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);

dimensionedVector gradPAMP2
(
    dimensionedVector::getOrDefault
    (
        "gradPAMP2",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);

dimensionedVector gradPAMP3
(
    dimensionedVector::getOrDefault
    (
        "gradPAMP3",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);

dimensionedVector gradPAMP4
(
    dimensionedVector::getOrDefault
    (
        "gradPAMP4",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);

dimensionedVector gradPAMP5
(
    dimensionedVector::getOrDefault
    (
        "gradPAMP5",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);

dimensionedVector gradPMEAN
(
    forceProperties.lookupOrDefault
    (
        "gradPMEAN",
        dimensionedVector("gradPMEAN",
                          dimensionSet(1, -2, -2, 0, 0, 0, 0),
                          vector(0, 0, 0) )
    )
);

dimensionedVector gradPOSC
(
    dimensionedVector::getOrDefault
    (
        "gradPOSC",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);

dimensionedVector gradConvAcc
(
    dimensionedVector::getOrDefault
    (
        "gradConvAcc",
        forceProperties,
        dimensionSet(1, -2, -2, 0, 0, 0, 0),
        vector(0, 0, 0)
    )
);
dimensionedScalar gradPOSCMin ("gradPOSCMin", gradPOSC.dimensions(), 1e-6);

dimensionedScalar oscpT
(
    forceProperties.lookupOrDefault
    (
        "oscpT",
        dimensionedScalar("oscpT",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0)
    )
);

// to avoid division by 0 in gradP impose 1 if 0 is specified
dimensionedScalar oscpTMin ("oscpTMin", oscpT.dimensions(), 1e-6);
oscpT=max(oscpTMin, oscpT);

dimensionedScalar Cvm
(
    forceProperties.lookupOrDefault
    (
        "Cvm",
        dimensionedScalar("Cvm",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0)
    )
);

dimensionedScalar Cl
(
    forceProperties.lookupOrDefault
    (
        "Cl",
        dimensionedScalar("Cl",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0)
    )
);

dimensionedScalar Ct
(
    forceProperties.lookupOrDefault
    (
        "Ct",
        dimensionedScalar("Ct",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0)
    )
);


dimensionedScalar initTheta1
(
    dimensionedScalar::getOrDefault
    (
        "initTheta1",
        forceProperties,
        dimless,
        0
    )
);

dimensionedScalar initTheta2
(
    dimensionedScalar::getOrDefault
    (
        "initTheta2",
        forceProperties,
        dimless,
        0
    )
);

dimensionedScalar initTheta3
(
    dimensionedScalar::getOrDefault
    (
        "initTheta3",
        forceProperties,
        dimless,
        0
    )
);

dimensionedScalar initTheta4
(
    dimensionedScalar::getOrDefault
    (
        "initTheta4",
        forceProperties,
        dimless,
        0
    )
);

dimensionedScalar initTheta5
(
    dimensionedScalar::getOrDefault
    (
        "initTheta5",
        forceProperties,
        dimless,
        0
    )
);
dimensionedScalar InitFreeze
(
    forceProperties.lookupOrDefault
    (
        "InitFreeze",
        dimensionedScalar("InitFreeze", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
    )
);
dimensionedScalar fixbeddepth
(
    forceProperties.lookupOrDefault
    (
        "fixbeddepth",
        dimensionedScalar("fixbeddepth", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0)
    )
);
dimensionedScalar ClipUa
(
    forceProperties.lookupOrDefault
    (
        "ClipUa",
        dimensionedScalar("ClipUa", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
    )
);

// Debuging info
Switch debugInfo
(
    forceProperties.lookupOrDefault("debugInfo", false)
);
// Sponge Layer
Switch spongeLayer
(
    forceProperties.lookupOrDefault("spongeLayer", false)
);
dimensionedScalar xSmin
(
    forceProperties.lookupOrDefault
    (
        "xSmin",
        dimensionedScalar("xSmin", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0)
    )
);
dimensionedScalar xSmax
(
    forceProperties.lookupOrDefault
    (
        "xSmax",
        dimensionedScalar("xSmax", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0)
    )
);

Switch writeTau
(
    forceProperties.lookupOrDefault<Switch>("writeTau", false)
);

Switch writeMomentumBudget
(
    forceProperties.getOrDefault<Switch>("writeMomentumBudget", false)
);



Info<< "Initializing with specified pressure gradient:"
<< gradPOSC.value() << endl << endl;

// open files to write gradP info
OFstream gradPOSCFile
(
    "gradPOSC.txt"
);

// Initialize variables
volScalarField Unity
(
    IOobject
    (
        "Unity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0), 1.0)
);

Info<< "Calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh
    ),
    fvc::interpolate(alpha)*phia + fvc::interpolate(beta)*phib
);

//#include "UaUbLimiter.H"
Ua.correctBoundaryConditions();

Info<< "Reading field p_rbgh\n" << endl;
volScalarField p_rbgh
(
    IOobject
    (
        "p_rbgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field pa\n" << endl;
volScalarField pa
(
    IOobject
    (
        "pa",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
);

Info<< "Reading field CohesionDistrb\n" << endl;

volScalarField CohesionDistrb
(
    IOobject
    (
        "CohesionDistrb",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volTensorField gradUaT
(
    fvc::grad(Ua)().T()
);

volTensorField gradUbT
(
    fvc::grad(Ub)().T()
);

volTensorField Taub
(
    IOobject
    (
        "Taub",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhob*(scalar(1.0)-alpha)*nub*dev(fvc::grad(Ub))
);

volTensorField Taua
(
    IOobject
    (
        "Taua",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    rhoa*alpha*nua*dev(fvc::grad(Ua))
);
//creating alphaPlastic field
Info<< "Computing field alphaPlastic\n" << endl;
volScalarField alphaPlastic
(
     IOobject
     (
         "alphaPlastic",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::NO_WRITE
     ),
     mesh,
     dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField delta
(
    IOobject
    (
        "delta",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
);


Info<< "Calculating field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha*Ua + beta*Ub
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    alpha*rhoa + beta*rhob
);

Info<< "Calculating field DDtUa and DDtUb\n" << endl;

volVectorField DDtUa
(
    fvc::ddt(Ua)
  + fvc::div(phia, Ua)
  - fvc::div(phia)*Ua
);

volVectorField DDtUb
(
    fvc::ddt(Ub)
  + fvc::div(phib, Ub)
  - fvc::div(phib)*Ub
);




Info<< "Calculating field g.h\n" << endl;

Info<< "\nReading hRef" << endl;
uniformDimensionedScalarField hRef
(
	IOobject
	(
		"hRef",
		runTime.constant(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	),
	dimensionedScalar(dimLength, Zero)
);
dimensionedScalar ghRef
(
	mag(gravity.value()) > SMALL
  ? gravity & (cmptMag(gravity.value())/mag(gravity.value()))*hRef
  : dimensionedScalar("ghRef", gravity.dimensions()*dimLength, 0)
);

volScalarField gh("gh", (g & mesh.C()) - ghRef);
surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);



// Read interfacial properties

IOdictionary interfacialProperties
(
    IOobject
    (
        "interfacialProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

autoPtr<dragModel> draga = dragModel::New
(
    interfacialProperties,
    phasea,
    phaseb
);

autoPtr<dragModel> dragb = dragModel::New
(
    interfacialProperties,
    phaseb,
    phasea
);

word dragPhase("blended");
if (interfacialProperties.found("dragPhase"))
{
    dragPhase = word(interfacialProperties.lookup("dragPhase"));

    bool validDrag =
        dragPhase == "a" || dragPhase == "b" || dragPhase == "blended";

    if (not validDrag)
    {
        FatalErrorIn(args.executable())
            << "invalid dragPhase " << dragPhase
            << exit(FatalError);
    }
}


Info << "dragPhase is " << dragPhase << endl;

// Compute fall velocity
#include "calculateUfall.H"
dimensionedScalar Ufall
(
   "Ufall",
   dimensionSet(0, 1, -1, 0, 0, 0, 0),
   max(UfallField).value()
);

// Initalize particle pressure model
autoPtr<ppModel> pp = ppModel::New
(
    ppProperties,
    phasea,
    phaseb
);

// Initalize Kinetic theory
kineticTheoryModel kineticTheory
(
    phasea,
    Ub,
    draga
);

// Initialize granular rheology
granularRheologyModel granularRheology
(
    phasea,
    phaseb,
    pa,
    Dsmall
);

// Declare variables
surfaceScalarField rUaAf
(
    IOobject
    (
        "rUaAf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 0.0)
);

volScalarField lambdaUa
(
    IOobject
    (
        "lambdaUa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
);

volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    draga->K(mag(Ub-Ua))
);

volVectorField gradPs
(
    IOobject
    (
        "gradPs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector::zero)
);

volVectorField gradPf
(
    IOobject
    (
        "gradPf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector::zero)
);
volVectorField susa
(
    IOobject
    (
        "susa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector::zero)
);
volVectorField susb
(
    IOobject
    (
        "susb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector::zero)
);
volScalarField pff
(
    IOobject
    (
        "pff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
);
volScalarField pS
(
    IOobject
    (
        "pS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pff + pa
);

Info<< "Calculating field muEff\n" << endl;
volScalarField muEff
(
    IOobject
    (
        "muEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1 , -1 ,-1 ,0 , 0), 0.0)
);
Info<< "Calculating field muFra\n" << endl;
volScalarField muFra
(
    IOobject
    (
        "muFra",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1 , -1 ,-1 ,0 , 0), 0.0)
);

volScalarField SUStilde
(
    IOobject
    (
        "SUStilde",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0), 1.0)
);

volVectorField SolidPressureContr
(
    IOobject
    (
        "SolidPressureContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(p_rbgh)
);

volVectorField divTauaContr
(
    IOobject
    (
        "divTauaContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);


volVectorField ViscStressContr
(
    IOobject
    (
        "ViscStressContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField FrictContr
(
    IOobject
    (
        "FrictContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField phiRaContr
(
    IOobject
    (
        "phiRaContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField phiRa2Contr
(
    IOobject
    (
        "phiRa2Contr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField divRcaContr
(
    IOobject
    (
        "divRcaContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField RcaContr
(
    IOobject
    (
        "RcaContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField BulkContr
(
    IOobject
    (
        "BulkContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField gravityContr
(
    IOobject
    (
        "gravityContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField DragContr
(
    IOobject
    (
        "DragContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField SuspContr
(
    IOobject
    (
        "SuspContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField ExtContr
(
    IOobject
    (
        "ExtContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    p_rbgh + rhob*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rbgh,
    mesh.solutionDict().subDict("PIMPLE"),
    pRefCell,
    pRefValue
);
if (p_rbgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rbgh = p - rhob*gh;
}
Info<<"end of createFields"<<endl;


//- Overset specific

// Add solver-specific interpolations
{
    wordHashSet& nonInt =
    const_cast<wordHashSet&>(Stencil::New(mesh).nonInterpolatedFields());
    nonInt.insert("surfaceIntegrate(phia)");
    nonInt.insert("surfaceIntegrate(phib)");
    nonInt.insert("cellMask");
    nonInt.insert("cellDisplacement");
    nonInt.insert("interpolatedCells");
    nonInt.insert("cellInterpolationWeight");

}
// Mask field for zeroing out contributions on hole cells
#include "createCellMask.H"

// Create bool field with interpolated cells
#include "createInterpolatedCells.H"
