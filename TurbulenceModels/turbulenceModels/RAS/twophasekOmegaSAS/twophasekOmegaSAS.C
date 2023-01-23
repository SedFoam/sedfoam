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

\*---------------------------------------------------------------------------*/

#include "twophasekOmegaSAS.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void twophasekOmegaSAS<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = k_/
    (
        max(omega_, Clim_*sqrt((2.0*(magSqr(symm(fvc::grad(this->U_)))))/Cmu_))
    );

    this->nut_.min(nutMax_);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> twophasekOmegaSAS<BasicTurbulenceModel>::Qsas
(
    const volScalarField::Internal& S2,
    const dimensionedScalar& alphaOmega,
    const dimensionedScalar& betaOmega
) const
{
        volScalarField::Internal L
    (
        sqrt(this->k_())/(pow025(this->Cmu_)*this->omega_())
    );

    volScalarField::Internal Lvk
    (
	max
        (
            kappa_*sqrt(S2)
           /(
                mag(fvc::laplacian(this->U_))()()
              + dimensionedScalar
                (
                    "ROOTVSMALL",
                    dimensionSet(0, -1, -1, 0, 0),
                    ROOTVSMALL
                )
            ),
            Cs_*sqrt(kappa_*zeta2_/((betaOmega_
            /this->Cmu_) - alphaOmega_))*this->delta()()
	)
    );

    return fvm::Su
    (
        this->alpha_()*this->rho_()
       *min
        (
            max
            (
                zeta2_*kappa_*S2*sqr(L/Lvk)
              - (2*C_/sigmaPhi_)*this->k_()
               *max
                (
                    magSqr(fvc::grad(this->omega_)()())/sqr(this->omega_()),
                    magSqr(fvc::grad(this->k_)()())/sqr(this->k_())
                ),
                dimensionedScalar(dimensionSet(0, 0, -2, 0, 0), Zero)
             ),
             // Limit SAS production of omega for numerical stability
             // particularly during start-up
             this->omega_()/(0.1*this->omega_.time().deltaT())
        ),
        this->omega_
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
twophasekOmegaSAS<BasicTurbulenceModel>::twophasekOmegaSAS
(
    const alphaField& beta,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& betaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        beta,
        rho,
        U,
        betaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    popeCorrection_
    (
        Switch::getOrAddToDict
        (
            "popeCorrection",
            this->coeffDict_,
            true
        )
    ),
    writeTke_
    (
        Switch::getOrAddToDict
        (
            "writeTke",
            this->coeffDict_,
            false
        )
    ),
    C3om_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C3om",
            this->coeffDict_,
            0.35
        )
    ),
    C4om_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C4om",
            this->coeffDict_,
            1.0
        )
    ),
    KE2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "KE2",
            this->coeffDict_,
            1.0
        )
    ),
    KE4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "KE4",
            this->coeffDict_,
            1.0
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    betaOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaOmega",
            this->coeffDict_,
            0.072
        )
    ),
    nutMax_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "nutMax",
            this->coeffDict_,
            0
        )
    ),
    Clim_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Clim",
            this->coeffDict_,
            0.0
        )
    ),
    sigmad_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmad",
            this->coeffDict_,
            0.0
        )
    ),
    alphaKOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaKOmega",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.52
        )
    ),
    alphaOmegaOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmegaOmega",
            this->coeffDict_,
            0.5
        )
    ),
    tmfexp_(U.db().lookupObject<volScalarField> ("tmfexp")),
    ESD3_(U.db().lookupObject<volScalarField> ("ESD3")),
    ESD4_(U.db().lookupObject<volScalarField> ("ESD4")),
    ESD5_(U.db().lookupObject<volScalarField> ("ESD5")),
    ESD_(U.db().lookupObject<volScalarField> ("ESD")),
    
    Cs_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.11
        )
    ),
    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    zeta2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "zeta2",
            this->coeffDict_,
            3.51
        )
    ),
    sigmaPhi_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaPhi",
            this->coeffDict_,
            2.0/3.0
        )
    ),
    C_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C",
            this->coeffDict_,
            2
        )
    ),

    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", betaRhoPhi.group()),
            *this,
            this->coeffDict_
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->correctNut();
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool twophasekOmegaSAS<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        betaOmega_.readIfPresent(this->coeffDict());
        alphaOmegaOmega_.readIfPresent(this->coeffDict());
        alphaKOmega_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
	Cs_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
        sigmaPhi_.readIfPresent(this->coeffDict());
        zeta2_.readIfPresent(this->coeffDict());
        C_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void twophasekOmegaSAS<BasicTurbulenceModel>::correct()
{
    if (not this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& beta = this->alpha_;
    const alphaField alpha = 1 - beta;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const surfaceScalarField& phi = this->phi_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    volSymmTensorField Sij(symm(fvc::grad(U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));

    volScalarField G
    (
        this->GName(),
        nut*2*magSqr(Sij)
    );

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volTensorField Omij(-skew(fvc::grad(U)));
    volVectorField Gradk(fvc::grad(k_));
    volVectorField Gradomega(fvc::grad(omega_));
    volScalarField alphadCheck_(Gradk & Gradomega);

    const volScalarField CDkOmega
    (
    sigmad_*pos(0.15-alpha)*pos(alphadCheck_)*(alphadCheck_)/omega_
    );


    volScalarField XsiOmega
    (
        IOobject
        (
            IOobject::groupName("XsiOmega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    if (popeCorrection_)
    {
        XsiOmega =
        (
            mag((Omij & Omij & Sij)/(pow((Cmu_*omega_), 3)))
        );
    }

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi, omega_)
      - fvm::Sp(fvc::div(phi), omega_)
      - fvm::laplacian(DomegaEff(), omega_)
      ==
      - fvm::SuSp (-alphaOmega_*G/k_, omega_)
      - fvm::Sp(ESD_, omega_)
      - fvm::Sp
        (
            betaOmega_*
            (
                (scalar(1.0)+scalar(85.0)*XsiOmega())
                /(scalar(1.0)+scalar(100.0)*XsiOmega())
            )*omega_(),
            omega_
        )
      + CDkOmega
      + ESD2()*fvm::Sp(C3om_*KE2_, omega_)
      + fvm::Sp((C4om_*KE4_*ESD5_*nut/k_), omega_)
      + Qsas(S2(), alphaOmega_, betaOmega_)
    );
    if (writeTke_)
    {
        #include "writeTKEBudget_kOmega.H"
    }
    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi, k_)
      - fvm::Sp(fvc::div(phi), k_)
      - fvm::laplacian(DkEff(), k_)
      ==
      - fvm::SuSp(-G/k_, k_)
      + fvm::Sp(-Cmu_*omega_, k_)
      + fvm::Sp(ESD_, k_)
      + fvm::Sp(KE4_*ESD4_*nut/k_, k_)
      + ESD2()*fvm::Sp(KE2_, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
