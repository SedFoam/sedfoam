/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::RASModels::twophasekOmega2006

Group
    grpRASTurbulence

Description
    Standard high Reynolds-number k-omega turbulence model for
    incompressible and compressible flows.

    References:
    \verbatim
        Wilcox, D. C. (1998).
        Turbulence modeling for CFD
        (Vol. 2, pp. 103-217). La Canada, CA: DCW industries.
    \endverbatim

    The default model coefficients are
    \verbatim
        twophasekOmega2006Coeffs
        {
            Cmu         0.09;  // Equivalent to betaStar
            alpha       0.52;
            betaOmega   0.072;
            alphak      0.5;
            alphaOmega  0.5;
        }
    \endverbatim

SourceFiles
    twophasekOmega2006.C

\*---------------------------------------------------------------------------*/

#ifndef twophasekOmega2006_H
#define twophasekOmega2006_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class twophasekOmega2006 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class twophasekOmega2006
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{

    protected:

    // Protected data

        // Dictionary

    IOdictionary twophaseRASProperties_;
    IOdictionary ppProperties_;
    // Model coefficients
    Switch popeCorrection_;
    Switch writeTke_;
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
    dimensionedScalar Clim_;
    dimensionedScalar alphaKOmega_;
    dimensionedScalar alphaOmega_;
    dimensionedScalar alphaOmegaOmega_;

    // Flow information
    // two phase information
    const volScalarField& tmfexp_;
    const volScalarField& ESD3_;
    const volScalarField& ESD4_;
    const volScalarField& ESD5_;
    const volScalarField& ESD_;
    dimensionedScalar alphaMax_;

    // Fields

    volScalarField k_;
    volScalarField omega_;

    // Protected Member Functions

    virtual void correctNut();


    public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("twophasekOmega2006");


    // Constructors

    //- Construct from components
    twophasekOmega2006
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName
    );


    //- Destructor
    #if (OFVERSION < 220 and OFVERSION >= 200) or OFVERSION == 00
        virtual ~twophasekOmega2006() = default;
    #else
        virtual ~twophasekOmega2006()
        {}
    #endif


    // Member Functions

    //- Read RASProperties dictionary
    virtual bool read();

    //- Return the effective diffusivity for k
    tmp<volScalarField> DkEff() const
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                "DkEff",
                alphaKOmega_*this->nut_ + this->nu()
            )
        );
    }

    //- Return the effective diffusivity for omega
    tmp<volScalarField> DomegaEff() const
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                "DomegaEff",
                alphaOmegaOmega_*this->nut_ + this->nu()
            )
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
                    this->mesh_.time().timeName(),
                    this->mesh_
                ),
                Cmu_*k_*omega_,
                omega_.boundaryField().types()
            )
        );
    }

    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "twophasekOmega2006.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
