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
    Foam::ChialvoSundaresanRadial

Description
    \verbatim
    "A modified kinetic theory for frictional granular flows in dense and 
    dilute regimes."
    Chialvo, S. and Sundaresan, S. (2013). Physics of Fluids, 25(7):070603.
    \endverbatim

SourceFiles
    ChialvoSundaresanRadial.C

\*---------------------------------------------------------------------------*/

#ifndef ChialvoSundaresanRadial_H
#define ChialvoSundaresanRadial_H

#include "radialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ChialvoSundaresanradial Declaration
\*---------------------------------------------------------------------------*/

class ChialvoSundaresanRadial
:
    public radialModel
{

public:

    //- Runtime type information
    TypeName("ChialvoSundaresan");


    // Constructors

        //- Construct from components
        ChialvoSundaresanRadial(const dictionary& dict);


    //- Destructor
    virtual ~ChialvoSundaresanRadial();


    // Member Functions

        tmp<volScalarField> g0
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax
        ) const;

        tmp<volScalarField> g0prime
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
