/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "Beetstra.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Beetstra, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        Beetstra,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Beetstra::Beetstra
(
    const dictionary& interfaceDict,
    const phaseModel& phasea,
    const phaseModel& phaseb
)
:
    dragModel(interfaceDict, phasea, phaseb)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Beetstra::~Beetstra()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::Beetstra::K
(
    const volScalarField& Ur
) const
{
    volScalarField Re(max(Ur*phasea_.d()/phaseb_.nu(), scalar(1.0e-3)));

    volScalarField F
    (
        10*alpha_/pow(1-alpha_, 2) +
        pow(1-alpha_, 2)*(1 + 1.5*pow(alpha_, 0.5)) +
        0.413*Re/(24*pow(1-alpha_, 2))*
        (pow(1-alpha_, -1) + 3*alpha_*(1-alpha_) + 8.4*pow(Re, -0.343))/
        (1+pow(10, 3*alpha_)*pow(Re, -(1+4*alpha_)/2))
    );
    volScalarField Cds
    (
        24/Re*pow(1-alpha_, -1)*F
    );

    return 0.75*Cds*phaseb_.rho()*Ur/phasea_.d();
}


// ************************************************************************* //
