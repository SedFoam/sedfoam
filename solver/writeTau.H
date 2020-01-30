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
if (runTime.outputTime())
{
    volTensorField Taub
    (
        IOobject
        (
            "Taub",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rhob*(scalar(1.0)-alpha)*nuEffb*dev(fvc::grad(Ub))
        );

    volTensorField Taua
    (
        IOobject
        (
            "Taua",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
        ),
        rhoa*(alpha*nuEffa+nuFra)*dev(fvc::grad(Ua))
    );

    runTime.write();
}
