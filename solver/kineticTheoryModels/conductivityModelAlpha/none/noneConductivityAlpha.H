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

Class
    Foam::noneViscosity

Description

SourceFiles
    noneViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef noneConductivityAlpha_H
#define noneConductivityAlpha_H

#include "conductivityModelAlpha.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class noneConductivityAlpha Declaration
\*---------------------------------------------------------------------------*/

class noneConductivityAlpha
:
    public conductivityModelAlpha
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noneConductivityAlpha(const dictionary& dict);


    //- Destructor
    virtual ~noneConductivityAlpha();


    // Member functions

        tmp<volScalarField> kappa
        (
            const volScalarField& alpha,
            const volScalarField& Theta,
            const volScalarField& g0,
            const volScalarField& g0Prime,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& da,
            const dimensionedScalar& e
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
