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

#include "DuanPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DuanPressure, 0);

    addToRunTimeSelectionTable
    (
        granularPressureModel,
        DuanPressure,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DuanPressure::DuanPressure(const dictionary& dict)
:
    granularPressureModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DuanPressure::~DuanPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::DuanPressure::granularPressureCoeff
(
    const volScalarField& alpha,
    const volScalarField& g0,
    //const volScalarField& Theta,
    const dimensionedScalar& rhoa,
    //const dimensionedScalar& da,
    const dimensionedScalar& e
) const
{
	//const dimensionedScalar kn = 5000;
	//const dimensionedScalar mu = 0.8;
	//const scalar sqrtPi = sqrt(constant::mathematical::pi);
	//const volScalarField tf = 1./(12./da*alpha*g0*sqrt(Theta)/sqrtPi);


    return rhoa*alpha*(1.0 + 2.0*(1.0 + e)*alpha*g0);
}


Foam::tmp<Foam::volScalarField> Foam::DuanPressure::granularPressureCoeffPrime
(
    const volScalarField& alpha,
    const volScalarField& g0,
    const volScalarField& g0prime,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& e
) const
{
    return rhoa*(1.0 + alpha*(1.0 + e)*(4.0*g0 + 2.0*g0prime*alpha));
}

// ************************************************************************* //
