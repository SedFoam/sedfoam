/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::incompressible::RASModels::twophasekEpsilon

Group
    grpIcoRASTurbulence

Description
    Standard k-epsilon turbulence model for incompressible flows.

    Reference:
    \verbatim
        Launder, B.E., and Spalding, D.B.,
        "Mathematical Models of Turbulence",
        Academic Press, 1972.

        Launder, B.E., and Spalding, D.B.,
        "Computational Methods for Turbulent Flows",
        Comp. Meth in Appl Mech & Eng'g, Vol 3, 1974, pp 269-289.
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        twophasekEpsilonCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            sigmaEps    1.3;
        }
    \endverbatim

SourceFiles
    twophasekEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef twophasekEpsilon_H
#define twophasekEpsilon_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class twophasekEpsilon Declaration
\*---------------------------------------------------------------------------*/

class twophasekEpsilon
:
    public RASModel
{

protected:

    // Protected data

          // Dictionary

            IOdictionary twophaseRASProperties_;
            IOdictionary ppProperties_;

        // Model coefficients

            dimensionedScalar C3ep_;
            dimensionedScalar C4ep_;
            dimensionedScalar KE2_;
            dimensionedScalar KE4_;
            dimensionedScalar omegaBC_;
            dimensionedScalar B_;

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar kSmall_;
            dimensionedScalar nutMax_;
            dimensionedScalar alphak_;
            dimensionedScalar alphaEps_;
            //dimensionedScalar sigmaEps_;

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
            volScalarField epsilon_;
            volScalarField nut_;


public:

    //- Runtime type information
    TypeName("twophasekEpsilon");

    // Constructors

        //- Construct from components
        twophasekEpsilon
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~twophasekEpsilon()
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
                new volScalarField("DkEff", nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/alphaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

          //- Return the fluid volume fraction
        tmp<volScalarField> beta() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("beta", 1.0 - alpha_)
            );
        }

            virtual tmp<volScalarField> omega() const
        {
	  //   return  epsilon_/(Cmu_*k_);
	    // /*
	  return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "omega",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    epsilon_/(Cmu_*k_),
                    epsilon_.boundaryField().types()
                )
            );
	 // */
        }

          virtual tmp<volScalarField> tmf() const
        {
	  return Foam::exp(tmfexp_*omega());
        }

        virtual tmp<volScalarField> ESD2() const
        {
            return ESD3_*(1.0-tmf());
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
