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
    Foam::MuIv

Description

SourceFiles
    MuIv.C

\*---------------------------------------------------------------------------*/

#ifndef MuIv_H
#define MuIv_H

#include "FrictionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace granularRheologyModels
{

/*---------------------------------------------------------------------------*\
                           Class MuIv Declaration
\*---------------------------------------------------------------------------*/

class MuIv:
public FrictionModel
{

    public:

    //- Runtime type information
    TypeName("MuIv");


    // Constructors

    //- Construct from components
    MuIv(const dictionary& dict);


    //- Destructor
    virtual ~MuIv();


    // Member functions

    tmp<volScalarField> muI
    (
        const dimensionedScalar& mus,
        const dimensionedScalar& mu2,
        const dimensionedScalar& I0,
        const volScalarField& pa,
        const dimensionedScalar& rhoa,
        const dimensionedScalar& da,
        const dimensionedScalar& rhob,
        const dimensionedScalar& nub,
        const volScalarField& magD,
        const dimensionedScalar& Dsmall
    ) const;

    tmp<volScalarField> I
    (
        const volScalarField& pa,
        const dimensionedScalar& rhoa,
        const dimensionedScalar& da,
        const dimensionedScalar& rhob,
        const dimensionedScalar& nub,
        const volScalarField& magD
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace granularRheologyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
