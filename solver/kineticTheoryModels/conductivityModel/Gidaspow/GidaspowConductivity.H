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
    Foam::GidaspowConductivity

Description
    \verbatim
    "A Bubbling Fluidization Model Using Kinetic Theory of Granular Flow",
    Ding, J., Gidaspow, D., AIChE Journal, 36(4):523–538
    \endverbatim

    This file is in the \dir kineticTheoryModels/conductivityModel/Gidaspow 
    directory

SourceFiles
    GidaspowConductivity.C

\*---------------------------------------------------------------------------*/

#ifndef GidaspowConductivity_H
#define GidaspowConductivity_H

#include "conductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class GidaspowConductivity Declaration
\*---------------------------------------------------------------------------*/

class GidaspowConductivity
:
    public conductivityModel
{

public:

        //- Runtime type information
        TypeName("Gidaspow");


    // Constructors

        //- Construct from components
        GidaspowConductivity(const dictionary& dict);


    //- Destructor
    virtual ~GidaspowConductivity();


    // Member Functions

        tmp<volScalarField> kappa
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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
