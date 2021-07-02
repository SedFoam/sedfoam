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
#include "readGravitationalAcceleration.H"

Info<< "\nStoring g as gravity\n";
uniformDimensionedVectorField gravity (g); // gMax already in use in CourantNo.H

// in newer OF versions, this should be expanded to also read g itself
IOdictionary gravitationalProperties
(
    IOobject
    (
        "gravitationalProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

const scalar gRampTime(gravitationalProperties.lookupOrDefault("gRampTime",
                                                                scalar(0)));

Info << "rampTime for gravity set to " << gRampTime << "\n";
