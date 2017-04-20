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

dimensionedScalar alphaSmall
(
    transportProperties.lookup("alphaSmall")
);

dimensionedScalar nuMax
(
    transportProperties.lookup("nuMax")
);
 
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

volVectorField& Ua = phasea->U();
surfaceScalarField& phia = phasea->phi();
const dimensionedScalar& rhoa = phasea->rho();
const dimensionedScalar& nua = phasea->nu();

volVectorField& Ub = phaseb->U();
surfaceScalarField& phib = phaseb->phi();
const dimensionedScalar& rhob = phaseb->rho();
const dimensionedScalar& nub = phaseb->nu();

Info<< "Reading field alpha\n" << endl;
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    scalar(1.0) - alpha
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
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
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
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

volVectorField DDtUa =
    fvc::ddt(Ua)
  + fvc::div(phia, Ua)
  - fvc::div(phia)*Ua;

volVectorField DDtUb =
    fvc::ddt(Ub)
  + fvc::div(phib, Ub)
  - fvc::div(phib)*Ub;


Info<< "Calculating field g.h\n" << endl;
volScalarField gh("gh", g & mesh.C());

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
    alpha,
    phasea,
    phaseb
);

autoPtr<dragModel> dragb = dragModel::New
(
    interfacialProperties,
    beta,
    phaseb,
    phasea
);

word dragPhase("blended");
if (interfacialProperties.found("dragPhase"))
{
    dragPhase = word(interfacialProperties.lookup("dragPhase"));

    bool validDrag =
        dragPhase == "a" || dragPhase == "b" || dragPhase == "blended";

    if (!validDrag)
    {
        FatalErrorIn(args.executable())
            << "invalid dragPhase " << dragPhase
            << exit(FatalError);
    }
}

Info << "dragPhase is " << dragPhase << endl;
kineticTheoryModel kineticTheory
(
    phasea,
    Ub,
    alpha,
    draga
);

granularRheologyModel granularRheology
(
    phasea,
    phaseb,
    alpha,
    draga,
    pa
);

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
    dimensionedScalar("zero", dimensionSet(1 , -1 ,-2 ,0 , 0), 0.0)
);

#include "createGradP.H"
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
