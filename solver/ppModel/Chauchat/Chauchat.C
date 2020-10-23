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

#include "Chauchat.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Chauchat, 0);

    addToRunTimeSelectionTable
    (
        ppModel,
        Chauchat,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Chauchat::Chauchat
(
    const dictionary& ppDict,
    const phaseModel& phasea,
    const phaseModel& phaseb
)
:
    ppModel(ppDict, phasea, phaseb)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Chauchat::~Chauchat()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::Chauchat::pff
(
    const volScalarField& alpha_,
    const volScalarField& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta0,
    const dimensionedScalar& eta1
) const
{
    scalar N = 2/(3-eta0.value());

    volScalarField arg1(1- max(alpha_ - alphaMinFriction, scalar(0))/alphaMax);

    return pos(alpha_-alphaMinFriction)*Fr*(pow(arg1, -N) - 1 );
}

// ************************************************************************* //
