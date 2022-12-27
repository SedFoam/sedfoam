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
    Foam::phaseModel

SourceFiles
    phaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseModel_H
#define phaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
{
    // Private data

        dictionary dict_;

        //- Name of phase
        word name_;

        //- Characteristic diameter of the solid phase
        dimensionedScalar d_;
        
        //- Characteristic diameter of large particles
        dimensionedScalar d_l_;

        //- shape factor for non-sperical particles
        dimensionedScalar sF_;
        
        //- 1st coefficient related to permeability for Engelund drag model
        dimensionedScalar aE_;
        
        //- 2nd coefficient related to permeability for Engelund drag model
        dimensionedScalar bE_;
        
        //- exponent of the hindrance function for drag coeeficient
        dimensionedScalar hExp_;

        //- kinematic viscosity
        dimensionedScalar nu_;

        //- density
        dimensionedScalar rho_;

        //- Velocity
        volVectorField U_;

        //- Concentration
        volScalarField alpha_;

        //- Fluxes
        autoPtr<surfaceScalarField> phiPtr_;


public:

    // Constructors

        phaseModel
        (
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const word& phaseName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<phaseModel> New
        (
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~phaseModel();


    // Member Functions

        const word& name() const
        {
            return name_;
        }

        const dimensionedScalar& d() const
        {
            return d_;
        }
        const dimensionedScalar& d_l() const
        {
            return d_l_;
        }
        const dimensionedScalar& sF() const
        {
            return sF_;
        }
        
        const dimensionedScalar& aE() const
        {
            return aE_;
        }
        
        const dimensionedScalar& bE() const
        {
            return bE_;
        }
        
        const dimensionedScalar& hExp() const
        {
            return hExp_;
        }

        const dimensionedScalar& nu() const
        {
            return nu_;
        }

        const dimensionedScalar& rho() const
        {
            return rho_;
        }

        const volVectorField& U() const
        {
            return U_;
        }

        volVectorField& U()
        {
            return U_;
        }

        const volScalarField& alpha() const
        {
            return alpha_;
        }

        volScalarField& alpha()
        {
            return alpha_;
        }

        const surfaceScalarField& phi() const
        {
            return phiPtr_();
        }

        surfaceScalarField& phi()
        {
            return phiPtr_();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
