/*---------------------------------------------------------------------------*\
Copyright (C) 2015 Cyrille Bonamy, Julien Chauchat, Tian-Jian Hsu
                   and contributors

License
    This file is part of SedFOAM.

    SedFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SedFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with SedFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "noneFriction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace granularRheologyModels
{
    defineTypeNameAndDebug(noneFriction, 0);
    addToRunTimeSelectionTable(FrictionModel, noneFriction, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::granularRheologyModels::noneFriction::noneFriction(const dictionary& dict)
:
    FrictionModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularRheologyModels::noneFriction::~noneFriction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::granularRheologyModels::noneFriction::muI
(
    const dimensionedScalar& mus,
    const dimensionedScalar& mu2,
    const dimensionedScalar& I0,
    const volScalarField& pa,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& da,
    const dimensionedScalar& rhob,
    const dimensionedScalar& nub,
    const volScalarField& magD,
    const dimensionedScalar& Dsmall
) const
{
    return scalar(0.0)*magD/(magD+Dsmall);
}

Foam::tmp<Foam::volScalarField> Foam::granularRheologyModels::noneFriction::I
(
    const volScalarField& pa,
    const dimensionedScalar& rhoa,
    const dimensionedScalar& da,
    const dimensionedScalar& rhob,
    const dimensionedScalar& nub,
    const volScalarField& magD
) const
{
    dimensionedScalar Dumsmall
    (
        "Dumsmall",
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        1e-4
    );
    return scalar(0.0)*magD/(magD+Dumsmall);
}
// ************************************************************************* //
