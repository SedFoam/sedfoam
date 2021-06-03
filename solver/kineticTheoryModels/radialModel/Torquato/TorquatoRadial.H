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
    Foam::TorquatoRadial

Description
    "Nearest-neighbor statistics for packings of hard spheres and disks."
    Torquato, S. (1995). Physical Review E, 51(4):3170–3182. 

    This file is in the \dir kineticTheoryModels/radialModel/Torquato directory

SourceFiles
    TorquatoRadial.C

\*---------------------------------------------------------------------------*/

#ifndef TorquatoRadial_H
#define TorquatoRadial_H

#include "radialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Torquatoradial Declaration
\*---------------------------------------------------------------------------*/

class TorquatoRadial
:
    public radialModel
{

public:

    //- Runtime type information
    TypeName("Torquato");


    // Constructors

        //- Construct from components
        TorquatoRadial(const dictionary& dict);


    //- Destructor
    virtual ~TorquatoRadial();


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
