/*---------------------------------------------------------------------------*\
Copyright (C) 2015 Cyrille Bonamy, Julien Chauchat, Tian-Jian Hsu
                   and contributors

License
    This file is part of SedFOAM.

    SedFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SedFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with SedFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::Chauchat

Description
    Chauchat, J., Guillou, S., Bang, D. P. V., and Nguyen, K. D. (2013).
    Modelling sedimentation–consolidation in the framework of a
    one-dimensional two-phase flow model. Journal of Hydraulic Research,
    51(3):293–305.

SourceFiles
    Chauchat.C

\*---------------------------------------------------------------------------*/

#ifndef Chauchat_H
#define Chauchat_H

#include "ppModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Chauchat Declaration
\*---------------------------------------------------------------------------*/

class Chauchat
:
    public ppModel
{

public:

    //- Runtime type information
    TypeName("Chauchat");


    // Constructors

        //- Construct from components
        Chauchat
        (
            const dictionary& ppDict,
            const phaseModel& phasea,
            const phaseModel& phaseb
        );


    //- Destructor
    virtual ~Chauchat();

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
