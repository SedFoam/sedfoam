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

#include "KriegerDoughertyViscosity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace granularRheologyModels
    {
        defineTypeNameAndDebug(KriegerDoughertyViscosity, 0);
        addToRunTimeSelectionTable
        (
            FluidViscosityModel, KriegerDoughertyViscosity, dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::granularRheologyModels::KriegerDoughertyViscosity::
                              KriegerDoughertyViscosity(const dictionary& dict):
    FluidViscosityModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularRheologyModels::KriegerDoughertyViscosity::
                             ~KriegerDoughertyViscosity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::granularRheologyModels::
                                      KriegerDoughertyViscosity::nuvb
(
    const volScalarField& alpha,
    const dimensionedScalar& nub,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Alphasmall,
    const dimensionedScalar& n
) const
{
    return nub*pow(1.0 - min(alpha/alphaMax, 0.99), -n)/(1.0-alpha);
}

// ************************************************************************* //
