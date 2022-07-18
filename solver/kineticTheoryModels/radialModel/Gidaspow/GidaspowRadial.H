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
    Foam::GidaspowRadial

Description
    \verbatim
    "A Bubbling Fluidization Model Using Kinetic Theory of Granular Flow",
    Ding, J., Gidaspow, D., AIChE Journal, 36(4):523–538
    \endverbatim

    This file is in the \dir kineticTheoryModels/radialModel/Gidaspow directory

SourceFiles
    GidaspowRadial.C

\*---------------------------------------------------------------------------*/

#ifndef GidaspowRadial_H
#define GidaspowRadial_H

#include "radialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Gidaspowradial Declaration
\*---------------------------------------------------------------------------*/

class GidaspowRadial
:
    public radialModel
{

public:

    //- Runtime type information
    TypeName("Gidaspow");


    // Constructors

        //- Construct from components
        GidaspowRadial(const dictionary& dict);


    //- Destructor
    virtual ~GidaspowRadial();


    // Member Functions

        tmp<volScalarField> g0
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax,
            const dimensionedScalar& muPart
        ) const;

        tmp<volScalarField> g0prime
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax,
            const dimensionedScalar& muPart
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
