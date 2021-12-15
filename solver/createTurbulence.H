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
#include "createTwoPhaseTurbulence.H"

singlePhaseTransportModel phaseaTransport(Ua, phia);

volScalarField alphaEff = alpha+10*alphaSmall;
surfaceScalarField alphacf("alphacf", fvc::interpolate(alphaEff));
surfaceScalarField alphaPhia("alphaPhia", alphacf*phia);

autoPtr<sedIncompressibleTurbulenceModel<singlePhaseTransportModel>>
turbulencea
(
    sedIncompressibleTurbulenceModel<singlePhaseTransportModel>::New
    (
        alphaEff,
        Ua,
        alphaPhia,
        phia,
        phaseaTransport,
        "turbulenceProperties.phasea"
    )
);

singlePhaseTransportModel phasebTransport(Ub, phib);

surfaceScalarField betacf("betacf", fvc::interpolate(beta));
surfaceScalarField betaPhib("betaPhib", betacf*phib);

autoPtr<sedIncompressibleTurbulenceModel<singlePhaseTransportModel>>
turbulenceb
(
    sedIncompressibleTurbulenceModel<singlePhaseTransportModel>::New
    (
        beta,
        Ub,
        betaPhib,
        phib,
        phasebTransport,
        "turbulenceProperties.phaseb"
    )
);

IOdictionary turbulencePropertiesa
(
    IOobject
    (
        "turbulenceProperties.phasea",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

IOdictionary turbulencePropertiesb
(
    IOobject
    (
        "turbulenceProperties.phaseb",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
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

volScalarField spherSigmaSGSb
(
    IOobject
    (
        "sphereSigmaSGSb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0), 0.0)
);

volScalarField spherSigmaSGSa
(
    IOobject
    (
        "spherSigmaSGSa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0), 0.0)
);
