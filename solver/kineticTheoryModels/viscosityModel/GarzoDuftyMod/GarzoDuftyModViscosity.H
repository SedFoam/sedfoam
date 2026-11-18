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
    Foam::GarzoDuftyModViscosity

Description
    Garzo and Dufty (1999) conductivity model modified empirically to account
    for saltation.
    \verbatim
    "A modified kinetic theory for frictional-collisional bedload transport
    valid from dense to dilute regime."
    Chassagne, R., Chauchat, J. and Bonamy, C. (2021).
    arXiv, physics.flu-dyn, https://arxiv.org/abs/2110.07340.
    \endverbatim
    \endverbatim

    This file is in the \dir kineticTheoryModels/viscosityModel/GarzoDuftyMod
    subdirectory inside the subdirectory
    \dir kineticTheoryModels/viscosityModel inside the directory
    \dir kineticTheoryModels

SourceFiles
    GarzoDuftyModViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef GarzoDuftyModViscosity_H
#define GarzoDuftyModViscosity_H

#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class GarzoDuftyModViscosity Declaration
\*---------------------------------------------------------------------------*/

class GarzoDuftyModViscosity
:
    public viscosityModel
{

public:

    //- Runtime type information
    TypeName("GarzoDuftyMod");


    // Constructors

        //- Construct from components
        GarzoDuftyModViscosity(const dictionary& dict);


    //- Destructor
    virtual ~GarzoDuftyModViscosity();


    // Member functions

        tmp<volScalarField> mua
        (
            const volScalarField& alpha,
            const volScalarField& Theta,
            const volScalarField& g0,
            const volScalarField& musalt,
            const volScalarField& K,
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
