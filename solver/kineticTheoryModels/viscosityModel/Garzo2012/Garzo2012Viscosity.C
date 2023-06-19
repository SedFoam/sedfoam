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

#include "Garzo2012Viscosity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(Garzo2012Viscosity, 0);
    addToRunTimeSelectionTable(viscosityModel, Garzo2012Viscosity, dictionary);
} // End namespace kineticTheoryModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::Garzo2012Viscosity::Garzo2012Viscosity
(
    const dictionary& dict
)
:
    viscosityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::Garzo2012Viscosity::~Garzo2012Viscosity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::Garzo2012Viscosity::mua
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
    const scalar pi = constant::mathematical::pi;

    //Impact of the deviation of the velocity distribution function
    //to the Maxwellian form (probably negligible)
    const dimensionedScalar mu4_0 = (7./2+pow(e, 2));
    const dimensionedScalar mu4_1 = (3./32*(69+10*pow(e, 2))+2./(1-e));
    const dimensionedScalar a2 = - (mu4_0-5)/(mu4_1-5*19./16);

    //Impact of fluid on the kinetic contribution
    const volScalarField gamma = (1-alpha)*K/rhoa;

    //Kinetic viscosity
    const volScalarField muk = alpha*pow(Theta, 0.5) *
            (1-2./5*(1+e)*(1-3*e)*alpha*g0) /
            (
             96/(5*sqrtPi)*(1-0.25*pow((1-e), 2)*(1+7./16*a2)- 5./24*
             (1-pow(e, 2))*(1+3*a2/16))*g0*alpha*pow(Theta, 0.5) +
             (1-alpha)*da*K/rhoa
            );
    //Contact viscosity
    const volScalarField muc = muk*(4./5*(1+e)*alpha*g0);
    //Bulk viscosity
    const volScalarField mub = 4./(5*sqrtPi)*(1+e)*(1-a2/16)*pow(alpha, 2)*g0;

    //Total viscosity accounting for saltation
    const volScalarField muTot = muk + muc + mub;

    return rhoa*da*sqrt(Theta)*muTot;
}

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::Garzo2012Viscosity::lambda
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
     1152./(45*pi)*(1+e)*pow(alpha, 2)*g0
    );
}

// ************************************************************************* //
