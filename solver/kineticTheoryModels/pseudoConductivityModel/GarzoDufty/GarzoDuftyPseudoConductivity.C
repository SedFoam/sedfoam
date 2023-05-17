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

#include "GarzoDuftyPseudoConductivity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GarzoDuftyPseudoConductivity, 0);

    addToRunTimeSelectionTable
    (
        pseudoConductivityModel,
        GarzoDuftyPseudoConductivity,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GarzoDuftyPseudoConductivity::GarzoDuftyPseudoConductivity
(
const dictionary& dict
)
:
    pseudoConductivityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GarzoDuftyPseudoConductivity::~GarzoDuftyPseudoConductivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GarzoDuftyPseudoConductivity::kappaAlpha
(
    const volScalarField& alpha,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& g0prime,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& da,
    const dimensionedScalar& e
) const
{
    const scalar sqrtPi = sqrt(constant::mathematical::pi);

    //Kinetic conductivity
    const volScalarField kappak = 125*sqrtPi/64 *
            (
             (1+3./5*pow(1+e, 2)*(2*e-1)*alpha*g0)/((1-7./16*(1-e))*(1+e)*g0)\
             *(1-pow(e, 2))*(g0+alpha*g0prime)
             - 6./25*alpha*(g0+alpha/2*g0prime)*e*(1-pow(e, 2))
            )/
            ((1+3./16*(1-e))*(1+e)*g0);
    //Contact conductivity
    const volScalarField kappac = 6./5*alpha*g0*(1+e)*kappak;

    //Total conductivity
    const volScalarField kappaTot = kappak+kappac;

    return rhoa*da*sqrt(Theta)*Theta*kappaTot;
}


// ************************************************************************* //
