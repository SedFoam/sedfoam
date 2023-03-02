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
    Foam::GarzoDuftyViscosity

Description
    \verbatim
    "Dense fluid transport for inelastic hard spheres."
    Garzó, V. and Dufty, J. W. (1999).
    Physical Review. E, Statistical Physics, Plasmas, Fluids, and Related
    Interdisciplinary Topics, 59(5 Pt B):5895–5911
    \endverbatim

    This file is in the \dir kineticTheoryModels/viscosityModel/GarzoDufty
    subdirectory inside the subdirectory
    \dir kineticTheoryModels/viscosityModel inside the directory
    \dir kineticTheoryModels

SourceFiles
    GarzoDuftyViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef GarzoDuftyViscosity_H
#define GarzoDuftyViscosity_H

#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class GarzoDuftyViscosity Declaration
\*---------------------------------------------------------------------------*/

class GarzoDuftyViscosity
:
    public viscosityModel
{

public:

    //- Runtime type information
    TypeName("GarzoDufty");


    // Constructors

        //- Construct from components
        GarzoDuftyViscosity(const dictionary& dict);


    //- Destructor
    virtual ~GarzoDuftyViscosity();


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
