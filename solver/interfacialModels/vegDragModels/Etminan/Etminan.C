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

#include "Etminan.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Etminan, 0);

    addToRunTimeSelectionTable
    (
        vegDragModel,
        Etminan,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Etminan::Etminan
(
    const dictionary& interfaceDict,
    const phaseModel& phasea,
    const phaseModel& phaseb
)
:
    vegDragModel(interfaceDict, phasea, phaseb)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Etminan::~Etminan()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::Etminan::K
(
    const volScalarField& Ua
) const
{

    const scalar Pi = constant::mathematical::pi;
    volScalarField Uc = (1-beta_)/(1-pow(2*beta_/Pi, 0.5))*Ua;
    volScalarField Rec
    (
         max(Uc*phaseb_.d()/phasea_.nu(), scalar(1.0e-9))
    );

    volScalarField Cd
    (
     1+10*pow(Rec, -2./3)
    );

    return (2.*Cd*phasea_.rho()/(Pi*alpha_*phaseb_.d())*
            pow((1-beta_)/(1-pow(2*beta_/Pi, 0.5)), 2)*Ua);
}


// ************************************************************************* //
