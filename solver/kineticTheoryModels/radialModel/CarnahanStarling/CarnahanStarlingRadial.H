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
    Foam::CarnahanStarlingRadial

Description
    \verbatim
    "Equation of State for Nonattracting Rigid Spheres."
    Carnahan, N. F. and Starling, K. E. (1969).
    The Journal of Chemical Physics, 51(2):635–636.
    \endverbatim

    This file is in the \dir kineticTheoryModels/radialModel/CarnahanStarling
    subdirectory inside the subdirectory
    \dir kineticTheoryModels/radialModel inside the directory
    \dir kineticTheoryModels

SourceFiles
    CarnahanStarlingRadial.C

\*---------------------------------------------------------------------------*/

#ifndef CarnahanStarlingRadial_H
#define CarnahanStarlingRadial_H

#include "radialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class CarnahanStarlingradial Declaration
\*---------------------------------------------------------------------------*/

class CarnahanStarlingRadial
:
    public radialModel
{

public:

    //- Runtime type information
    TypeName("CarnahanStarling");


    // Constructors

        //- Construct from components
        CarnahanStarlingRadial(const dictionary& dict);


    //- Destructor
    virtual ~CarnahanStarlingRadial();


    // Member Functions


        tmp<volScalarField> g0
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax
        ) const;

        tmp<volScalarField> g0prime
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
