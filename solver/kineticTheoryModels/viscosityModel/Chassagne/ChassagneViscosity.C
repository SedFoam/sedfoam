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

#include "ChassagneViscosity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(ChassagneViscosity, 0);
    addToRunTimeSelectionTable(viscosityModel, ChassagneViscosity, dictionary);
} // End namespace kineticTheoryModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::ChassagneViscosity::ChassagneViscosity
(
    const dictionary& dict
)
:
    viscosityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::ChassagneViscosity::~ChassagneViscosity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::kineticTheoryModels::ChassagneViscosity::mua
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
    const scalar pi = constant::mathematical::pi;

    return rhoa*da*sqrt(Theta)*5*sqrtPi/96*
    (
     //Kinetic viscosity
     (48/(5*sqrtPi)*alpha-2./5*(1+e)*(1-3*e)*alpha*g0)/((1-0.25*pow((1-e),2)-5./24*(1-pow(e,2)))*g0)+
     //Contact viscosity
     (48/(5*sqrtPi)*alpha-2./5*(1+e)*(1-3*e)*alpha*g0)/((1-0.25*pow((1-e),2)-5./24*(1-pow(e,2)))*g0)*(4./5*(1+e)*alpha*g0) +
     //Bulk viscosity
     384./(25*pi)*(1+e)*pow(alpha,2)*g0
    );
}


// ************************************************************************* //
