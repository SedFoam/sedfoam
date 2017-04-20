/*---------------------------------------------------------------------------*\
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
    Foam::incompressible::RASModels::twophasekOmega

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
        twophasekOmegaCoeffs
        {
            Cmu         0.09;  // Equivalent to betaStar
            alpha       0.52;
            betaOmega   0.072;
            alphak      0.5;
            alphaOmega  0.5;
        }
    \endverbatim

SourceFiles
    twophasekOmega.C

\*---------------------------------------------------------------------------*/

#ifndef twophasekOmega_H
#define twophasekOmega_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class twophasekOmega Declaration
\*---------------------------------------------------------------------------*/

class twophasekOmega
:
    public RASModel
{

//FOAMEXT COMM protected:

    // Protected data

        // Dictionary

            IOdictionary twophaseRASProperties_;
            IOdictionary ppProperties_;

        // Model coefficients

            dimensionedScalar C3om_;
            dimensionedScalar C4om_;
            dimensionedScalar KE2_;
            dimensionedScalar KE4_;
            dimensionedScalar omegaBC_;
            dimensionedScalar B_;

            dimensionedScalar Cmu_;
            dimensionedScalar betaOmega_;
            dimensionedScalar kSmall_;
            dimensionedScalar nutMax_;
//            dimensionedScalar alphaK_;
            dimensionedScalar alphaKOmega_;
            dimensionedScalar alphaOmega_;
            dimensionedScalar alphaOmegaOmega_;




        // flow information
            const volScalarField& alpha_;
        // two phase information
            const volScalarField& tmfexp_;
            const volScalarField& ESD3_;
            const volScalarField& ESD4_;
            const volScalarField& ESD5_;
            const volScalarField& ESD_;
    
            scalar alphaMax_;
            // Fields

            volScalarField k_;
            volScalarField omega_;
            volScalarField nut_;

public:

    //- Runtime type information
    TypeName("twophasekOmega");

    // Constructors

        //- Construct from components
        twophasekOmega
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~twophasekOmega()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphaKOmega_*nut_ + nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff", alphaOmegaOmega_*nut_ + nu())
            );
        }

        //- Return the fluid volume fraction
        tmp<volScalarField> beta() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("beta", 1.0 - alpha_)
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        virtual tmp<volScalarField> tmf() const
        {
            return Foam::exp(tmfexp_*omega_);
        }

        virtual tmp<volScalarField> ESD2() const
        {
            return ESD3_*(1.0-tmf());
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    Cmu_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
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
