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
    Foam::RASModels::twophasekOmegaSAS

Group
    grpRASTurbulence

Description
    Standard high Reynolds-number k-omega SAS turbulence model for
    incompressible and compressible flows.

    References:
    \verbatim
        k-omega SST-SAS original model:
            Egorov, Y. (2008).
            Development and Application of SST-SAS Turbulence Model
            in the DESIDER Project
            Advances in Hybrid RANS-LES Modelling

        k-omega original model:
            Wilcox, D. C. (1998).
            Turbulence modeling for CFD
            (Vol. 2, pp. 103-217). La Canada, CA: DCW industries.

        Two-phase adaptation:
            Amoudry, L. O. (2014).
            Extension of k-omega turbulence closure to two-phase sediment
            transport modelling: Application to oscillatory sheet flows.
            Advances in Water Ressources, 72, 110-121
    \endverbatim

    The default model coefficients are
    \verbatim
        twophasekOmegaSASCoeffs
        {
        // Default kOmega Wilcox (1998) coefficients
            Cmu              0.09;
            alphaOmega       0.52;
            betaOmega        0.072;
            C3om             0.35;
            C4om             1.0;
            KE2              1.0;
            KE4              1.0;
            alphaKomega      0.5;
            alphaOmegaOmega  0.5;
            Clim             0.0;
            sigmad           0.0
            popeCorrection   false;
        }
    \endverbatim

    For the revisited version of the model from Wilcox (2006),
    model coefficients should be modified following:
    \verbatim
        twophasekOmegaSASCoeffs
        {
        // Default kOmega Wilcox (2006) coefficients
            betaOmega        0.0708;
            alphaKomega      0.6;
            Clim             0.875;
            sigmad           0.125
            popeCorrection   True (only for 3D cases)

        // Default SAS coefficients
            Cs          0.11;
            kappa       0.41;
            zeta2       3.51;
            sigmaPhi    2.0/3.0;
            C           2;

            // Delta must be specified for SAS e.g.
            delta cubeRootVol;

            cubeRootVolCoeffs
            {}
        }
    \endverbatim

    This file is in the \dir turbulenceModels/RAS/twophasekOmegaSAS subdirectory
    inside the \dir turbulenceModels/RAS subdirectory
    inside the \dir turbulenceModels subdirectory

SourceFiles
    twophasekOmegaSAS.C

\*---------------------------------------------------------------------------*/

#ifndef twophasekOmegaSAS_H
#define twophasekOmegaSAS_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class twophasekOmegaSAS Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class twophasekOmegaSAS
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        twophasekOmegaSAS(const twophasekOmegaSAS&);
        void operator=(const twophasekOmegaSAS&);

protected:

    // Protected data

        // Model coefficients
            Switch popeCorrection_;
            Switch writeTke_;
            dimensionedScalar C3om_;
            dimensionedScalar C4om_;
            dimensionedScalar KE2_;
            dimensionedScalar KE4_;

            dimensionedScalar Cmu_;
            dimensionedScalar betaOmega_;
            dimensionedScalar nutMax_;
            dimensionedScalar Clim_;
            dimensionedScalar sigmad_;
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

          // Model constants

            dimensionedScalar Cs_;
            dimensionedScalar kappa_;
            dimensionedScalar zeta2_;
            dimensionedScalar sigmaPhi_;
            dimensionedScalar C_;

        // Fields
        volScalarField k_;
            volScalarField omega_;

            //- Run-time selectable delta model
                autoPtr<Foam::LESdelta> delta_;

    // Protected Member Functions

        //- SAS omega source
    virtual tmp<fvScalarMatrix> Qsas
        (
            const volScalarField::Internal& S2,
            const dimensionedScalar& alphaOmega,
            const dimensionedScalar& betaOmega
        ) const;

    virtual void correctNut();

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("twophasekOmegaSAS");


    // Constructors

        //- Construct from components
        twophasekOmegaSAS
        (
            const alphaField& beta,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& betaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~twophasekOmegaSAS()
    {}


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
            return Foam::exp(Cmu_*tmfexp_*omega_);
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

        //- Access function to filter width
        inline const volScalarField& delta() const
        {
            return *delta_;
        }

    //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "twophasekOmegaSAS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
