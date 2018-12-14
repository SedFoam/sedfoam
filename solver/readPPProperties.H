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

scalar alphaMax(readScalar(ppProperties.lookup("alphaMax")));
dimensionedScalar alphaMinFriction(ppProperties.lookup("alphaMinFriction"));
dimensionedScalar Fr_(ppProperties.lookup("Fr"));
dimensionedScalar eta0_
(
    ppProperties.lookupOrDefault
    (
       "eta0_",
       dimensionedScalar("eta0_", dimensionSet(0, 0, 0, 0, 0, 0, 0), 3 )
    )
);
dimensionedScalar eta1_
(
    ppProperties.lookupOrDefault
    (
       "eta1_",
       dimensionedScalar("eta1_", dimensionSet(0, 0, 0, 0, 0, 0, 0), 5 )
    )
);
Switch packingLimiter(ppProperties.lookup("packingLimiter"));
