/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "JohnsonJackson.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(JohnsonJackson, 0);

    addToRunTimeSelectionTable
    (
        ppModel,
        JohnsonJackson,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::JohnsonJackson::JohnsonJackson
(
    const dictionary& ppDict,
    const phaseModel& phasea,
    const phaseModel& phaseb
)
:
    ppModel(ppDict, phasea, phaseb)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::JohnsonJackson::~JohnsonJackson()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::JohnsonJackson::pff
(
    const volScalarField& alpha_,
    const volScalarField& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta0,
    const dimensionedScalar& eta1
) const
{
// limiting parameter : 1e-3 because :
// max(pff)=Fr * (0.625-0.53)^3/limiteur^5
// limiteur =(Fr (0.625-0.53)^3/max(pff))^(1/5)
// with max(pff)=0.6(rhop-rhof)max(hbed)
// limiteur =0.05
    return Fr*pow(max(alpha_ - alphaMinFriction, scalar(0)), eta0)
             /pow(max(alphaMax - alpha_, scalar(1.0e-3)), eta1);
}

// ************************************************************************* //
