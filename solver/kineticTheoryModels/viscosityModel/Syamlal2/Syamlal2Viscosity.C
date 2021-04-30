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

#include "Syamlal2Viscosity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(Syamlal2Viscosity, 0);
    addToRunTimeSelectionTable(viscosityModel, Syamlal2Viscosity, dictionary);
} // End namespace kineticTheoryModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::Syamlal2Viscosity::Syamlal2Viscosity
(
    const dictionary& dict
)
:
    viscosityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::Syamlal2Viscosity::~Syamlal2Viscosity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::kineticTheoryModels::Syamlal2Viscosity::mua
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

    return rhoa*da*sqrt(Theta)*
    (
        0.8*(4.0/5.0)*sqr(alpha)*g0*(1.0 + e)/sqrtPi
      + 0.8*(1.0/15.0)*sqrtPi*g0*(1.0 + e)*(3.0*e - 1.0)*sqr(alpha)/(3.0 - e)
      + (4.0/6.0)*alpha*sqrtPi/(3.0 - e)
    );
}


// ************************************************************************* //
