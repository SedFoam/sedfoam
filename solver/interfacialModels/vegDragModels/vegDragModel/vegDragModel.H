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
    Foam::vegDragModel

Description

SourceFiles
    vegDragModel.C
    newVegDragModel.C

\*---------------------------------------------------------------------------*/

#ifndef vegDragModel_H
#define vegDragModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "phaseModel.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class vegDragModel Declaration
\*---------------------------------------------------------------------------*/

class vegDragModel
{
protected:

    // Protected data

    const dictionary& interfaceDict_;
    const volScalarField& alpha_;
    const volScalarField& beta_;
    const phaseModel& phasea_;
    const phaseModel& phaseb_;


public:

    //- Runtime type information
    TypeName("vegDragModel");


    // Declare runtime construction

    declareRunTimeSelectionTable
    (
        autoPtr,
        vegDragModel,
        dictionary,
        (
            const dictionary& interfaceDict,
            const phaseModel& phasea,
            const phaseModel& phaseb
        ),
        (interfaceDict, phasea, phaseb)
    );


    // Constructors

    vegDragModel
    (
        const dictionary& interfaceDict,
        const phaseModel& phasea,
        const phaseModel& phaseb
    );


    //- Destructor
    virtual ~vegDragModel();


    // Selectors

    static autoPtr<vegDragModel> New
    (
        const dictionary& interfaceDict,
        const phaseModel& phasea,
        const phaseModel& phaseb
    );


    // Member Functions

    //- the dragfunction K used in the momentum eq.
    //    ddt(alpha.rhoa.Ua) + ... = ... -alpha.beta.K.Ua
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
    // for numerical reasons alpha and beta has been
    // extracted from the dragFunction K,
    // so you MUST divide K by alpha.beta when implementing the vegetation
    // drag function
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
    virtual tmp<volScalarField> K(const volScalarField& U) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
