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
    Foam::Etminan

Description
    Etminan, V., R. J. Lowe, and M. Ghisalberti (2017),
    A new model for predicting the drag exerted by vegetation canopies,
    Water Resour. Res., 53, 31793196, doi:10.1002/ 2016WR020090.

SourceFiles
    Etminan.C

\*---------------------------------------------------------------------------*/

#ifndef Etminan_H
#define Etminan_H

#include "vegDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DallaValle Declaration
\*---------------------------------------------------------------------------*/

class Etminan
:
    public vegDragModel
{

public:

    //- Runtime type information
    TypeName("Etminan");


    // Constructors

        //- Construct from components
        Etminan
        (
            const dictionary& interfaceDict,
            const phaseModel& phasea,
            const phaseModel& phaseb
        );


    //- Destructor
    virtual ~Etminan();


    // Member Functions

        tmp<volScalarField> K(const volScalarField& Ua) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
