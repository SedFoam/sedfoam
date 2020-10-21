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
    Foam::noneDilatancy

Description

SourceFiles
    noneDilatancy.C

\*---------------------------------------------------------------------------*/

#ifndef noneDilatancy_H
#define noneDilatancy_H

#include "DilatancyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace granularRheologyModels
{

/*---------------------------------------------------------------------------*\
                           Class noneDilatancy Declaration
\*---------------------------------------------------------------------------*/

class noneDilatancy
:
    public DilatancyModel
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noneDilatancy(const dictionary& dict);


    //- Destructor
    virtual ~noneDilatancy();


    // Member functions

    tmp<volScalarField> delta
    (
    const dimensionedScalar& K_dila,
    const dimensionedScalar& alpha_c,
    const volScalarField& alpha,
    const volScalarField& magD,
    const dimensionedScalar& da,
    const dimensionedScalar& rhob,
    const dimensionedScalar& nub,
    const volScalarField& p_p_total,
    const dimensionedScalar& PaMin
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace granularRheologyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
