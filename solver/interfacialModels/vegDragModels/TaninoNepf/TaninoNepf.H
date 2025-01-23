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
    Foam::TaninoNepf

Description
    Tanino, Y., & Nepf, H. M. (2008a). Laboratory investigation of mean drag
    in a random array of rigid, emergent cylinders. Journal of Hydraulic
    Engineering, 134(1), 3441.
    https://doi.org/10.1061/(ASCE)0733-9429(2008)134:1(34)

SourceFiles
    TaninoNepf.C

\*---------------------------------------------------------------------------*/

#ifndef TaninoNepf_H
#define TaninoNepf_H

#include "vegDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DallaValle Declaration
\*---------------------------------------------------------------------------*/

class TaninoNepf
:
    public vegDragModel
{

public:

    //- Runtime type information
    TypeName("TaninoNepf");


    // Constructors

        //- Construct from components
        TaninoNepf
        (
            const dictionary& interfaceDict,
            const phaseModel& phasea,
            const phaseModel& phaseb
        );


    //- Destructor
    virtual ~TaninoNepf();


    // Member Functions

        tmp<volScalarField> K(const volScalarField& Ua) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
