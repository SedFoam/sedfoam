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
    Foam::LunPressure

Description
    "Kinetic theories for granular flow: Inelastic particles in Couette flow 
    and slightly inelastic particles in a general flowfield." 
    Lun, C. K. K., Savage, S. B., Jeffrey, D. J., and Chepurniy, N. (1984). 
    Journal of Fluid Mechanics, 140:223–256.

SourceFiles
    LunPressure.C

\*---------------------------------------------------------------------------*/

#ifndef LunPressure_H
#define LunPressure_H

#include "granularPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LunPressure Declaration
\*---------------------------------------------------------------------------*/

class LunPressure
:
    public granularPressureModel
{

public:

    //- Runtime type information
    TypeName("Lun");


    // Constructors

        //- Construct from components
        LunPressure(const dictionary& dict);


    //- Destructor
    virtual ~LunPressure();


    // Member Functions

        tmp<volScalarField> granularPressureCoeff
        (
            const volScalarField& alpha,
            const volScalarField& g0,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& e
        ) const;

        tmp<volScalarField> granularPressureCoeffPrime
        (
            const volScalarField& alpha,
            const volScalarField& g0,
            const volScalarField& g0prime,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& e
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
