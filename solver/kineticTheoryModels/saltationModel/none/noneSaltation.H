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
    Foam::noneSaltation

Description
    This file is in the \dir kineticTheoryModels/saltationModel/none directory

SourceFiles
    noneSaltation.C

\*---------------------------------------------------------------------------*/

#ifndef noneSaltation_H
#define noneSaltation_H

#include "saltationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class noneSaltation Declaration
\*---------------------------------------------------------------------------*/

class noneSaltation
:
    public saltationModel
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

    //- Construct from components
    noneSaltation(const dictionary& dict);


    //- Destructor
    virtual ~noneSaltation();


    // Member functions

        tmp<volScalarField> musalt
        (
            const volScalarField& alpha,
            const volScalarField& Theta,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& da,
            const volScalarField& K
        ) const;

        tmp<volScalarField> kappasalt
        (
            const volScalarField& alpha,
            const volScalarField& Theta,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& da,
            const volScalarField& K
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
