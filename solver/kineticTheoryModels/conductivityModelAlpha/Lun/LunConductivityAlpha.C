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

#include "LunConductivityAlpha.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LunConductivityAlpha, 0);

    addToRunTimeSelectionTable
    (
        conductivityModelAlpha,
        LunConductivityAlpha,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LunConductivityAlpha::LunConductivityAlpha(const dictionary& dict)
:
    conductivityModelAlpha(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LunConductivityAlpha::~LunConductivityAlpha()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::LunConductivityAlpha::kappa
(
    const volScalarField& alpha,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& g0Prime,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& da,
    const dimensionedScalar& e
) const
{
    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    const dimensionedScalar r = 1./2*(1+e);

    return rhoa*da*sqrt(Theta)*Theta*
    (
     15./4*(2*r-1)*(r-1)*sqrtPi/(41-33*r)*(1./(alpha*g0)+12./5*r)*
     (2*alpha*g0 + alpha*alpha*g0Prime)
    );
}


// ************************************************************************* //
