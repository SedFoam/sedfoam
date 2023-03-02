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
    Foam::Garzo2012Conductivity

Description
    \verbatim
    "Enskog kinetic theory for monodisperse gas–solid flows."
    Garzó, V., Tenneti, S., Subramaniam, S. & Hrenya, C. M. (2012).
    Journal of Fluid Mechanics 712, 129–168.
    \endverbatim

    This file is in the \dir kineticTheoryModels/conductivityModel/Garzo2012
    directory

SourceFiles
    Garzo2012Conductivity.C

\*---------------------------------------------------------------------------*/

#ifndef Garzo2012Conductivity_H
#define Garzo2012Conductivity_H

#include "conductivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Garzo2012Conductivity Declaration
\*---------------------------------------------------------------------------*/

class Garzo2012Conductivity
:
    public conductivityModel
{

public:

    //- Runtime type information
    TypeName("Garzo2012");


    // Constructors

        //- Construct from components
        Garzo2012Conductivity(const dictionary& dict);


    //- Destructor
    virtual ~Garzo2012Conductivity();


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
