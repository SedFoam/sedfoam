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
    Foam::HrenyaSinclairConductivity

Description
    This file is in the \dir
    kineticTheoryModels/conductivityModel/HrenyaSinclair directory

SourceFiles
    HrenyaSinclairConductivity.C

\*---------------------------------------------------------------------------*/

#ifndef HrenyaSinclairConductivity_H
#define HrenyaSinclairConductivity_H

#include "conductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class HrenyaSinclairConductivity Declaration
\*---------------------------------------------------------------------------*/

class HrenyaSinclairConductivity
:
    public conductivityModel
{

        dictionary coeffsDict_;

        //- characteristic length of geometry
        dimensionedScalar L_;

public:

        //- Runtime type information
        TypeName("HrenyaSinclair");


    // Constructors

        //- Construct from components
        HrenyaSinclairConductivity(const dictionary& dict);


    //- Destructor
    virtual ~HrenyaSinclairConductivity();


    // Member Functions

        tmp<volScalarField> kappa
        (
            const volScalarField& alpha,
            const volScalarField& Theta,
            const volScalarField& g0,
            const volScalarField& kappasalt,
            const volScalarField& K,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& da,
            const dimensionedScalar& e
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
