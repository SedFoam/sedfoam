/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::incompressible::RASModels::twophaseMixingLength

Group
    grpIcoRASTurbulence

Description
    Standard high Reynolds-number k-omega turbulence model for
    incompressible flows.

    References:
    \verbatim
        "Turbulence Modeling for CFD"
        D. C. Wilcox,
        DCW Industries, Inc., La Canada,
        California, 1988.

        See also:
        http://www.cfd-online.com/Wiki/Wilcox's_k-omega_model
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        twophaseMixingLengthCoeffs
        {
            Cmu         0.09;  // Equivalent to betaStar
            alpha       0.52;
            betaOmega   0.072;
            alphak      0.5;
            alphaOmega  0.5;
        }
    \endverbatim

SourceFiles
    twophaseMixingLength.C

\*---------------------------------------------------------------------------*/

#ifndef twophaseMixingLength_H
#define twophaseMixingLength_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class twophaseMixingLength Declaration
\*---------------------------------------------------------------------------*/

class twophaseMixingLength
:
    public RASModel
{

//FOAMEXT COMM protected:

    // Protected data

        // Dictionary

            IOdictionary twophaseRASProperties_;

        // Model coefficients
            dimensionedScalar alphaMaxLM_;
            dimensionedScalar kappaLM_;

        // flow information
            volScalarField alpha_;

        // Fields
            volScalarField k_;
            volScalarField epsilon_;
            volScalarField nut_;

public:

    //- Runtime type information
    TypeName("twophaseMixingLength");

    // Constructors

        //- Construct from components
        twophaseMixingLength
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~twophaseMixingLength()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
