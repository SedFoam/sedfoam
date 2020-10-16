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

dimensionedScalar SUSread
(
    twophaseRASProperties.lookupOrDefault
    (
        "SUS",
        dimensionedScalar("SUS", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
    )
);
dimensionedScalar SUS_I
(
    twophaseRASProperties.lookupOrDefault
    (
        "SUS_I",
        dimensionedScalar("SUS_I", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1)
    )
);

dimensionedScalar SUS_A
(
    twophaseRASProperties.lookupOrDefault
    (
        "SUS_A",
        dimensionedScalar("SUS_A", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
    )
);

volScalarField SUS
(
    IOobject
    (
        "SUS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    SUSread
);

Switch SUSlocal
(
    twophaseRASProperties.lookupOrDefault("SUSlocal", false)
);

dimensionedScalar Usmall
(
    transportProperties.lookupOrDefault
    (
        "Usmall",
        dimensionedScalar("Usmall", dimensionSet(0, 1, -1, 0, 0, 0, 0), 1e-6)
    )
);

dimensionedScalar KE1
(
    twophaseRASProperties.lookupOrDefault
    (
        "KE1",
        dimensionedScalar("KE1", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
    )
);
dimensionedScalar Tpsmall
(
    twophaseRASProperties.lookupOrDefault
    (
        "Tpsmall",
        dimensionedScalar("Tpsmall", dimensionSet(0, 0, 1, 0, 0, 0, 0), 1e-6)
    )
);
dimensionedScalar KE3
(
    twophaseRASProperties.lookupOrDefault
    (
        "KE3",
        dimensionedScalar("KE3", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
    )
);
dimensionedScalar Cmu
(
    twophaseRASProperties.lookupOrDefault
    (
        "Cmu",
        dimensionedScalar("Cmu", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.09)
    )
);
dimensionedScalar B
(
    twophaseRASProperties.lookupOrDefault
    (
        "B",
        dimensionedScalar("B", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.25)
    )
);
volScalarField tmfexp
(
    IOobject
    (
        "tmfexp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -B*rhoa*6.0*Cmu/(max(beta*K, Tpsmall))
);

volScalarField rtemp = (Ub-Ua) & gravity;
volVectorField Urtemp = Ub-Ua;
if (Foam::mag(gravity).value() != 0)
{
	Urtemp = rtemp*gravity/magSqr(gravity);
}

volVectorField Urh = (Ub-Ua) - Urtemp;
volScalarField ESD1 = (K/rhob)*(fvc::grad(alpha)) & Urh;

volScalarField ESD3
(
    IOobject
    (
        "ESD3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -2.0*K*alpha/rhob
);

volScalarField ESD4
(
    IOobject
    (
        "ESD4",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (-(SUS/beta)*(rhoa/rhob-1.0)*gravity & fvc::grad(alpha)) + KE1*ESD1
);

volScalarField ESD5
(
    IOobject
    (
        "ESD5",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    max
    (
        -(SUS/beta)*(rhoa/rhob-1.0)*gravity & fvc::grad(alpha),
        scalar(0.0)*ESD1
    )
    + KE1*ESD1
);

volScalarField ESD
(
    IOobject
    (
        "ESD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    KE3*ESD3
);
