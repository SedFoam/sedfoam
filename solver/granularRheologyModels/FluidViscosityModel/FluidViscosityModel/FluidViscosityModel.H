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
    Foam::FluidViscosityModel

Description

SourceFiles
    FluidViscosityModel.C

\*---------------------------------------------------------------------------*/

#ifndef FluidViscosityModel_H
#define FluidViscosityModel_H

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
                           Class FluidViscosityModel Declaration
\*---------------------------------------------------------------------------*/

class FluidViscosityModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        FluidViscosityModel(const FluidViscosityModel&);

        //- Disallow default bitwise assignment
        void operator=(const FluidViscosityModel&);


protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("FluidViscosityModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        FluidViscosityModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        FluidViscosityModel(const dictionary& dict);


    // Selectors

        static autoPtr<FluidViscosityModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~FluidViscosityModel();


    // Member Functions

        virtual tmp<volScalarField> nuvb
        (
	 const volScalarField& alpha,
	 const dimensionedScalar& nub,
	 const dimensionedScalar& alphaMax,
	 const dimensionedScalar& Alphasmall,
	 const dimensionedScalar& n
        )const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace granularRheologyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
