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

#include "Coulomb.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace granularRheologyModels
{
    defineTypeNameAndDebug(Coulomb, 0);
    addToRunTimeSelectionTable(FrictionModel, Coulomb, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::granularRheologyModels::Coulomb::Coulomb(const dictionary& dict)
:
    FrictionModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularRheologyModels::Coulomb::~Coulomb()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::granularRheologyModels::Coulomb::muI
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
    return mus*(magD+Dsmall)/(magD+Dsmall);
// mus is multiplied by a volScalarField to be consistent with the general
// definition of muI. (magD+Dsmall)/(magD+Dsmall) is chosen because it is always
// egal to one independently from the value of magD even if it is stricly zero.
}

Foam::tmp<Foam::volScalarField> Foam::granularRheologyModels::Coulomb::I
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
    return 0*magD/(magD+Dumsmall);
}

// ************************************************************************* //
