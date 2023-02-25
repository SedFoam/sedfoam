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
    Foam::LunSavageRadial

Description
    "The effects of an impact velocity dependent coefficient of restitution
    on stresses developed by sheared granular materials."
    Lun, C. K. K. and Savage, S. B. 1986. Acta Mech., 63, 15–44.

    This file is in the \dir kineticTheoryModels/radialModel/LunSavage
    subdirectory inside the subdirectory
    \dir kineticTheoryModels/radialModel inside the directory
    \dir kineticTheoryModels

SourceFiles
    LunSavageRadial.C

\*---------------------------------------------------------------------------*/

#ifndef LunSavageRadial_H
#define LunSavageRadial_H

#include "radialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LunSavageradial Declaration
\*---------------------------------------------------------------------------*/

class LunSavageRadial
:
    public radialModel
{

public:

    //- Runtime type information
    TypeName("LunSavage");


    // Constructors

        //- Construct from components
        LunSavageRadial(const dictionary& dict);


    //- Destructor
    virtual ~LunSavageRadial();


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
