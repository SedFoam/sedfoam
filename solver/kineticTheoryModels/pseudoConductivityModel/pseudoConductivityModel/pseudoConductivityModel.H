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
    Foam::pseudoConductivityModel

Description
    \verbatim
    "Dense fluid transport for inelastic hard spheres."
    Garzó, V. and Dufty, J. W. (1999).
    Physical Review. E, Statistical Physics, Plasmas, Fluids, and Related
    Interdisciplinary Topics, 59(5 Pt B):5895–5911
    \endverbatim

SourceFiles
    pseudoConductivityModel.C

\*---------------------------------------------------------------------------*/

#ifndef pseudoConductivityModel_H
#define pseudoConductivityModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pseudoConductivityModel Declaration
\*---------------------------------------------------------------------------*/

class pseudoConductivityModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        pseudoConductivityModel(const pseudoConductivityModel&);

        //- Disallow default bitwise assignment
        void operator=(const pseudoConductivityModel&);


protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("pseudoConductivityModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        pseudoConductivityModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        pseudoConductivityModel(const dictionary& dict);


    // Selectors

        static autoPtr<pseudoConductivityModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~pseudoConductivityModel();


    // Member Functions

        virtual tmp<volScalarField> kappaAlpha
        (
            const volScalarField& alpha,
            const volScalarField& Theta,
            const volScalarField& g0,
            const volScalarField& g0prime,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& da,
            const dimensionedScalar& e
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
