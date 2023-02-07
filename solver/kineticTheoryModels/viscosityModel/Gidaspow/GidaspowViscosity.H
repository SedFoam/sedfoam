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
    Foam::GidaspowViscosity

Description
    \verbatim
    "A Bubbling Fluidization Model Using Kinetic Theory of Granular Flow",
    Ding, J., Gidaspow, D., AIChE Journal, 36(4):523–538
    \endverbatim

    This file is in the \dir kineticTheoryModels/viscosityModel/Gidaspow
    subdirectory inside the subdirectory
    \dir kineticTheoryModels/viscosityModel inside the directory
    \dir kineticTheoryModels

SourceFiles
    GidaspowViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef GidaspowViscosity_H
#define GidaspowViscosity_H

#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class GidaspowViscosity Declaration
\*---------------------------------------------------------------------------*/

class GidaspowViscosity
:
    public viscosityModel
{

public:

    //- Runtime type information
    TypeName("Gidaspow");

    // Constructors

        //- Construct from components
        GidaspowViscosity(const dictionary& dict);


    //- Destructor
    virtual ~GidaspowViscosity();


    // Member functions

        tmp<volScalarField> mua
        (
            const volScalarField& alpha,
            const volScalarField& Theta,
            const volScalarField& g0,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& da,
            const dimensionedScalar& e
        ) const;

        tmp<volScalarField> lambda
        (
            const volScalarField& alpha,
            const volScalarField& Theta,
            const volScalarField& g0,
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
