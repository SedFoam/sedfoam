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
    Foam::MerckelbachKranenburg

Description
    Chauchat, J., Guillou, S., Bang, D. P. V., and Nguyen, K. D. (2013).
    Modelling sedimentation–consolidation in the framework of a
    one-dimensional two-phase flow model. Journal of Hydraulic Research,
    51(3):293–305.

SourceFiles
    MerckelbachKranenburg.C

\*---------------------------------------------------------------------------*/

#ifndef MerckelbachKranenburg_H
#define MerckelbachKranenburg_H

#include "ppModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class MerckelbachKranenburg Declaration
\*---------------------------------------------------------------------------*/

class MerckelbachKranenburg
:
    public ppModel
{

public:

    //- Runtime type information
    TypeName("MerckelbachKranenburg");


    // Constructors

        //- Construct from components
        MerckelbachKranenburg
        (
            const dictionary& ppDict,
            const phaseModel& phasea,
            const phaseModel& phaseb
        );


    //- Destructor
    virtual ~MerckelbachKranenburg();

    // Member Functions

    tmp<volScalarField> pff
    (
        const volScalarField& alpha_,
        const volScalarField& alphaMinFriction,
        const dimensionedScalar& alphaMax,
        const dimensionedScalar& Fr,
        const dimensionedScalar& eta0,
        const dimensionedScalar& eta1
      ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
