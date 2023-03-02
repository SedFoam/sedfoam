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
    Foam::saltationModel

Description
    This file is in the \dir
    kineticTheoryModels/saltationModel/saltationModel directory

SourceFiles
    saltationModel.C

\*---------------------------------------------------------------------------*/

#ifndef saltationModel_H
#define saltationModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class saltationModel Declaration
\*---------------------------------------------------------------------------*/

class saltationModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        saltationModel(const saltationModel&);

        //- Disallow default bitwise assignment
        void operator=(const saltationModel&);


protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("saltationModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        saltationModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        saltationModel(const dictionary& dict);


    // Selectors

        static autoPtr<saltationModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~saltationModel();


    // Member Functions

        virtual tmp<volScalarField> musalt
        (
            const volScalarField& alpha,
            const volScalarField& Theta,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& da,
            const volScalarField& K
        ) const = 0;

	virtual tmp<volScalarField> kappasalt
        (
            const volScalarField& alpha,
            const volScalarField& Theta,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& da,
            const volScalarField& K
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
