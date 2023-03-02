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

#include "HrenyaSinclairViscosity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(HrenyaSinclairViscosity, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        HrenyaSinclairViscosity,
        dictionary
    );
} // End namespace kineticTheoryModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::HrenyaSinclairViscosity::HrenyaSinclairViscosity
(
    const dictionary& dict
)
:
    viscosityModel(dict),
#if (defined(OPENFOAM) && (OPENFOAM >= 1812))
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs")),
    L_("L", dimLength, coeffsDict_)
#else
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    L_(coeffsDict_.get<scalar>("L"))
#endif
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::HrenyaSinclairViscosity::~HrenyaSinclairViscosity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::HrenyaSinclairViscosity::mua
(
    const volScalarField& alpha,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& musalt,
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
        (4.0/5.0)*sqr(alpha)*g0*(1.0 + e)/sqrtPi
      + (1.0/15.0)*sqrtPi*g0*(1.0 + e)*(3.0*e - 1)*sqr(alpha)/(3.0-e)
      + (1.0/6.0)*sqrtPi*alpha*(0.5*lamda + 0.25*(3.0*e - 1.0))
       /(0.5*(3.0 - e)*lamda)
      + (10/96.0)*sqrtPi/((1.0 + e)*0.5*(3.0 - e)*g0*lamda)
    );
}

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::HrenyaSinclairViscosity::lambda
(
    const volScalarField& alpha,
    const volScalarField& Theta,
    const volScalarField& g0,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& da,
    const dimensionedScalar& e
) const
{
    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    // Lun and Savage bulk viscosity (to improve)
    return rhoa*da*sqrt(Theta)*
    (
        (4.0/3.0)*sqr(alpha)*g0*(1.0+e)/sqrtPi
    );
}
// ************************************************************************* //
