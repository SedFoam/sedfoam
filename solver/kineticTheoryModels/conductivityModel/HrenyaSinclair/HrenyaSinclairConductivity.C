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

#include "HrenyaSinclairConductivity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HrenyaSinclairConductivity, 0);

    addToRunTimeSelectionTable
    (
        conductivityModel,
        HrenyaSinclairConductivity,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HrenyaSinclairConductivity::HrenyaSinclairConductivity
(
    const dictionary& dict
)
:
    conductivityModel(dict),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    L_
    (
        dimensionedScalar::getOrDefault
        (
            "L",
            coeffsDict_,
            dimensionSet(0, 1, 0, 0, 0, 0, 0),
            0.0005
         )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HrenyaSinclairConductivity::~HrenyaSinclairConductivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::HrenyaSinclairConductivity::kappa
(
    const volScalarField& alpha,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& kappasalt,
    const volScalarField& K,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& da,
    const dimensionedScalar& e
) const
{
    const scalar sqrtPi = sqrt(constant::mathematical::pi);

    volScalarField lamda
    (
        scalar(1) + da/(6.0*sqrt(2.0)*(alpha + scalar(1.0e-5)))/L_
    );

    return rhoa*da*sqrt(Theta)*
    (
        2.0*sqr(alpha)*g0*(1.0 + e)/sqrtPi
      + (9.0/8.0)*sqrtPi*0.25*sqr(1.0 + e)*(2.0*e - 1.0)*sqr(alpha)
       /(49.0/16.0 - 33.0*e/16.0)
      + (15.0/16.0)*sqrtPi*alpha*(0.5*sqr(e) + 0.25*e - 0.75 + lamda)
       /((49.0/16.0 - 33.0*e/16.0)*lamda)
      + (25.0/64.0)*sqrtPi
       /((1.0 + e)*(49.0/16.0 - 33.0*e/16.0)*lamda*g0)
    );
}


// ************************************************************************* //
