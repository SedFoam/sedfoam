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

#include "GarzoDuftyModConductivity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GarzoDuftyModConductivity, 0);

    addToRunTimeSelectionTable
    (
        conductivityModel,
        GarzoDuftyModConductivity,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GarzoDuftyModConductivity::
GarzoDuftyModConductivity(const dictionary& dict)
:
    conductivityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GarzoDuftyModConductivity::~GarzoDuftyModConductivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GarzoDuftyModConductivity::kappa
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

    //Kinetic conductivity
    const volScalarField kappak = 25*sqrtPi/64*
            (64/(25*sqrtPi)*alpha+3./5*pow(1+e, 2)*(2*e-1)*alpha*g0)/
            ((1-7/16*(1-e))*(1+e)*g0);
    //Contact conductivity
    const volScalarField kappac = kappak*6./5*(1+e)*alpha*g0;
    //Bulk conductivity
    const volScalarField kappab = 2/sqrtPi*(1+e)*pow(alpha, 2)*g0;

    //Total conductivity accounting for saltation
    const volScalarField kappaTot = kappak + kappac + kappab;

    return rhoa*da*sqrt(Theta)*kappaTot;
}


// ************************************************************************* //
