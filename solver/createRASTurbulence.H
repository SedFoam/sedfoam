/*---------------------------------------------------------------------------*\
Copyright (C) 2015 Cyrille Bonamy, Julien Chauchat, Tian-Jian Hsu
                   and contributors

License
    This file is part of SedFOAM.

    SedFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SedFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with SedFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/
#include "createTwoPhaseRASTurbulence.H"

singlePhaseTransportModel phasebTransport(Ub, phib);

surfaceScalarField alphacf("alphacf", fvc::interpolate(alpha));
surfaceScalarField alphaPhic("alphaPhic", alphacf*phib);

autoPtr<sedIncompressibleTurbulenceModel<singlePhaseTransportModel>>
turbulenceb
(
    sedIncompressibleTurbulenceModel<singlePhaseTransportModel>::New
    (
        alpha,
        Ub,
        alphaPhic,
        phib,
        phasebTransport
    )
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
        IOobject::AUTO_WRITE
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
        IOobject::AUTO_WRITE
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
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
);
