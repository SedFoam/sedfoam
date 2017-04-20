#include "createTwoPhaseRASTurbulence.H"

    singlePhaseTransportModel phasebTransport(Ub, phib);

    autoPtr<incompressible::turbulenceModel> turbulenceb
    (
        incompressible::turbulenceModel::New(Ub, phib, phasebTransport)
    );

Info<< "Calculating field nuEffa\n" << endl;
volScalarField nuEffa
(
    IOobject
    (
        "nuEffa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
);
Info<< "Calculating field nuFra\n" << endl;
volScalarField nuFra
(
    IOobject
    (
        "nuFra",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
);

Info<< "Calculating field nuEffb\n" << endl;
volScalarField nuEffb
(
    IOobject
    (
        "nuEffb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
);
