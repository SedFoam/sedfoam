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

#include "BoyerEtAlViscosity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace granularRheologyModels
{
    defineTypeNameAndDebug(BoyerEtAlViscosity, 0);
    addToRunTimeSelectionTable
    (
        FluidViscosityModel,
        BoyerEtAlViscosity,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::granularRheologyModels::BoyerEtAlViscosity::BoyerEtAlViscosity
(
    const dictionary& dict
):
    FluidViscosityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularRheologyModels::BoyerEtAlViscosity::~BoyerEtAlViscosity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::granularRheologyModels::
                                      BoyerEtAlViscosity::nuvb
(
    const volScalarField& alpha,
    const dimensionedScalar& nub,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Alphasmall,
    const dimensionedScalar& n
) const
{
    return nub*(1.0 + 2.5*alpha*pow
                                (
                                    1.0 - min(alpha/alphaMax, scalar(0.99)),
                                    -1
                                ))/(1.0-alpha);
}


// ************************************************************************* //
