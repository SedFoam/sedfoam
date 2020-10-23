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
    Foam::DilatancyModel

Description

SourceFiles
    DilatancyModel.C

\*---------------------------------------------------------------------------*/

#ifndef DilatancyModel_H
#define DilatancyModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace granularRheologyModels
{

/*---------------------------------------------------------------------------*\
                           Class DilatancyModel Declaration
\*---------------------------------------------------------------------------*/

class DilatancyModel
{
    // Private member functions

    //- Disallow default bitwise copy construct
    DilatancyModel(const DilatancyModel&);

    //- Disallow default bitwise assignment
    void operator=(const DilatancyModel&);


protected:

    // Protected data

    const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("DilatancyModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DilatancyModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

    // Construct from components
    DilatancyModel(const dictionary& dict);


    // Selectors

    static autoPtr<DilatancyModel> New
    (
        const dictionary& dict
    );


    //- Destructor
    virtual ~DilatancyModel();


    // Member Functions

    virtual tmp<volScalarField> delta
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
    ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace granularRheologyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
