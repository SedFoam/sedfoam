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

#include "Garzo2012Conductivity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Garzo2012Conductivity, 0);

    addToRunTimeSelectionTable
    (
        conductivityModel,
        Garzo2012Conductivity,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Garzo2012Conductivity::Garzo2012Conductivity(const dictionary& dict)
:
    conductivityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Garzo2012Conductivity::~Garzo2012Conductivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::Garzo2012Conductivity::kappa
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
    const scalar Pi = constant::mathematical::pi;

    //Impact of the deviation of the velocity distribution function
    //to the Maxwellian form (probably negligible)
    const volScalarField mu2_0 = pow(2*Pi, 0.5)*(1-pow(e, 2))*g0;
    const volScalarField mu2_1 = 3./16*mu2_0;
    const volScalarField mu4_0 = (7./2+pow(e, 2))*mu2_0;
    const volScalarField mu4_1 = (3./32*(69+10*pow(e, 2))+2./(1-e))*mu2_0;
    const volScalarField a2 = - (mu4_0-5*mu2_0)/(mu4_1-5*19./16*mu2_0);

    //Impact of fluid on the kinetic contribution
    const volScalarField gamma = (1-alpha)*K/rhoa;


    //Kinetic conductivity
    const volScalarField kappak = 25./64*sqrtPi*
            (1+2*a2 + 3./5*pow(1+e, 2)*(2*e-1+a2*(1+e)))/
            ((1-7./16*(1-e)+(827-459*e)/256*a2)*(1+e)*g0);
    //Contact conductivity
    const volScalarField kappac = kappak*6./5*(1+e)*alpha*g0;
    //Bulk conductivity
    const volScalarField kappab = 2/sqrtPi*(1+e)*(1+7/16*a2)*pow(alpha, 2)*g0;

    //Total conductivity accounting for saltation
    const volScalarField kappaTot = kappak + kappac + kappab;

    return rhoa*da*sqrt(Theta)*kappaTot;
}


// ************************************************************************* //
