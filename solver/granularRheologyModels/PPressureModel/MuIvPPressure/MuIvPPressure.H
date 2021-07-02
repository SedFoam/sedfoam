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
    Foam::MuIvPPressure

Description

    This file is in the \dir PPressureModel/MuIvPPressure subdirectory
    inside the \dir PPressureModel subdirectory

SourceFiles
    MuIvPPressure.C

\*---------------------------------------------------------------------------*/

#ifndef MuIvPPressure_H
#define MuIvPPressure_H

#include "PPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace granularRheologyModels
{

/*---------------------------------------------------------------------------*\
                           Class MuIvPPressure Declaration
\*---------------------------------------------------------------------------*/

class MuIvPPressure
:
    public PPressureModel
{

public:

    //- Runtime type information
    TypeName("MuIv");


    // Constructors

        //- Construct from components
        MuIvPPressure(const dictionary& dict);


    //- Destructor
    virtual ~MuIvPPressure();


    // Member functions
    tmp<volScalarField> pa
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
    ) const;

    tmp<volScalarField> alphaEq
    (
        const volScalarField& pa,
        const dimensionedScalar& Bphi,
        const dimensionedScalar& rhoa,
        const dimensionedScalar& da,
        const dimensionedScalar& rhob,
        const dimensionedScalar& nub,
        const volScalarField& magD,
        const dimensionedScalar& alphaMax
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace granularRheologyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
