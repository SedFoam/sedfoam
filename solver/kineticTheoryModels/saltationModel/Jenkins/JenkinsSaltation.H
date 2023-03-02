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
    Foam::JenkinsSaltation

Description
    \verbatim
    "Continuum model for steady, fully developed saltation above a horizontal
    particle bed"
    Jenkins, J. T., Cantat, I. & Valance, A. (2010).
    Physical Review E 82 (2), 020301.

    "Two-phase continuum theory for windblown sand."
    Jenkins, James T. & Valance, Alexandre (2018).
    Physical Review Fluids 3 (3), 034305.
    \endverbatim

    This file is in the \dir kineticTheoryModels/saltationModel/Jenkins
    directory

SourceFiles
    JenkinsSaltation.C

\*---------------------------------------------------------------------------*/

#ifndef JenkinsSaltation_H
#define JenkinsSaltation_H

#include "saltationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class JenkinsSaltation Declaration
\*---------------------------------------------------------------------------*/

class JenkinsSaltation
:
    public saltationModel
{

public:

    //- Runtime type information
    TypeName("Jenkins");


    // Constructors

        //- Construct from components
        JenkinsSaltation(const dictionary& dict);


    //- Destructor
    virtual ~JenkinsSaltation();


    // Member functions

	tmp<volScalarField> musalt
        (
            const volScalarField& alpha,
            const volScalarField& Theta,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& da,
            const volScalarField& K
        ) const;

        tmp<volScalarField> kappasalt
        (
            const volScalarField& alpha,
            const volScalarField& Theta,
            const dimensionedScalar& rhoa,
            const dimensionedScalar& da,
            const volScalarField& K
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
