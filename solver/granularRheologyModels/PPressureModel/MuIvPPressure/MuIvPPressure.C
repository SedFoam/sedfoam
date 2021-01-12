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

#include "MuIvPPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace granularRheologyModels
{
    defineTypeNameAndDebug(MuIvPPressure, 0);
    addToRunTimeSelectionTable(PPressureModel, MuIvPPressure, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::granularRheologyModels::MuIvPPressure::MuIvPPressure
(
    const dictionary& dict
)
:
    PPressureModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularRheologyModels::MuIvPPressure::~MuIvPPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::granularRheologyModels::MuIvPPressure::pa
(
    const volScalarField& pf,
    const dimensionedScalar& Bphi,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& da,
    const dimensionedScalar& rhob,
    const dimensionedScalar& nub,
    const volScalarField& magD,
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Alphasmall
) const
{
    // MuI viscous regime: pa = (b alpha / (alphaMax - alpha))^2 nub rhob magD
    return pow(Bphi*alpha / max(alphaMax-alpha, scalar(1e-3)), 2)*nub*rhob*magD;
}

Foam::tmp<Foam::volScalarField> Foam::granularRheologyModels::
MuIvPPressure::alphaEq
(
    const volScalarField& pa,
    const dimensionedScalar& Bphi,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& da,
    const dimensionedScalar& rhob,
    const dimensionedScalar& nub,
    const volScalarField& magD,
    const dimensionedScalar& alphaMax
) const

{
    //
    //  inertial regime: alphaMax / (1+Bphi I)

    return alphaMax/(1+Bphi*rhob*nub*magD/pa);
}
// ************************************************************************* //
