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

#include "ChialvoSundaresanRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ChialvoSundaresanRadial, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        ChialvoSundaresanRadial,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ChialvoSundaresanRadial::ChialvoSundaresanRadial(const dictionary& dict)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ChialvoSundaresanRadial::~ChialvoSundaresanRadial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::ChialvoSundaresanRadial::g0
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax
) const
{

    return neg(alpha-alphaMax)*min(1e4, (2-alpha)/(2*pow(1-alpha,3)) + 4.69*pow(alpha,2)/(alphaMax-alpha))
	    + pos(alpha-alphaMax)*1e4;
}


Foam::tmp<Foam::volScalarField> Foam::ChialvoSundaresanRadial::g0prime
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax
) const
{
    return 2./alphaMax*pow(mag(1.0 - alpha/alphaMax), -3.0);
}


// ************************************************************************* //
