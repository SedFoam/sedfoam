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

#include "noneSaltation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(noneSaltation, 0);
    addToRunTimeSelectionTable(saltationModel, noneSaltation, dictionary);
} // End namespace kineticTheoryModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::noneSaltation::noneSaltation(const dictionary& dict)
:
    saltationModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::noneSaltation::~noneSaltation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::kineticTheoryModels::noneSaltation::musalt
(
 	const volScalarField& alpha,
        const volScalarField& Theta,
        const dimensionedScalar& rhoa,
        const dimensionedScalar& da,
        const volScalarField& K
) const
{

        return dimensionedScalar
        (
            "1000000000000000",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            1000000000000000
        )*alpha/alpha;
}

Foam::tmp<Foam::volScalarField> Foam::kineticTheoryModels::noneSaltation::kappasalt
(
 	const volScalarField& alpha,
        const volScalarField& Theta,
        const dimensionedScalar& rhoa,
        const dimensionedScalar& da,
        const volScalarField& K
) const
{

        return dimensionedScalar
        (
            "1000000000000000",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            1000000000000000
        )*alpha/alpha;
}


// ************************************************************************* //
