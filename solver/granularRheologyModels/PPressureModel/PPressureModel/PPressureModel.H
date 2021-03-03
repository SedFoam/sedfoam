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
    Foam::PPressureModel

Description

SourceFiles
    PPressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef PPressureModel_H
#define PPressureModel_H

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
                           Class PPressureModel Declaration
\*---------------------------------------------------------------------------*/

class PPressureModel
{
    // Private member functions

    //- Disallow default bitwise copy construct
    PPressureModel(const PPressureModel&);

    //- Disallow default bitwise assignment
    void operator=(const PPressureModel&);


    protected:

    // Protected data

    const dictionary& dict_;


    public:

    //- Runtime type information
    TypeName("PPressureModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PPressureModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

    // Construct from components
    PPressureModel(const dictionary& dict);


    // Selectors

    static autoPtr<PPressureModel> New
    (
        const dictionary& dict
    );


    //- Destructor
    virtual ~PPressureModel();


    // Member Functions

    virtual tmp<volScalarField> pa
    (
        const volScalarField& pf,
        const dimensionedScalar& Bphi,
        const dimensionedScalar& rhoa,
        const dimensionedScalar& da,
        const dimensionedScalar& rhob,
        const dimensionedScalar& nub,
        const volScalarField& magD,
        const volScalarField& alpha,
        const dimensionedScalar& alphaMax,
        const dimensionedScalar& Alphasmall
    ) const = 0;

    virtual tmp<volScalarField> alphaEq
    (
        const volScalarField& pa,
        const dimensionedScalar& Bphi,
        const dimensionedScalar& rhoa,
        const dimensionedScalar& da,
        const dimensionedScalar& rhob,
        const dimensionedScalar& nub,
        const volScalarField& magD,
        const dimensionedScalar& alphaMax
    ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace granularRheologyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
