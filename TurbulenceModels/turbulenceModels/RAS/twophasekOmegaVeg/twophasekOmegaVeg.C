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

#include "twophasekOmegaVeg.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void twophasekOmegaVeg<BasicTurbulenceModel>::correctNut()
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
twophasekOmegaVeg<BasicTurbulenceModel>::twophasekOmegaVeg
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
    KE6_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "KE6",
            this->coeffDict_,
            0.2
        )
    ),
    KE7_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "KE7",
            this->coeffDict_,
            0.15
        )
    ),
    KE8_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "KE8",
            this->coeffDict_,
            1
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
    Clambda_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Clambda",
            this->coeffDict_,
            0.01
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
            1e-1
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
    ESD6_(U.db().lookupObject<volScalarField> ("ESD6")),
    ESD7_(U.db().lookupObject<volScalarField> ("ESD7")),
    ESD8_(U.db().lookupObject<volScalarField> ("ESD8")),
    ESD_(U.db().lookupObject<volScalarField> ("ESD")),

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
    ks_
    (
        IOobject
        (
            IOobject::groupName("ks", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    kw_
    (
        IOobject
        (
            IOobject::groupName("kw", U.group()),
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
    ),
    omegas_
    (
        IOobject
        (
            IOobject::groupName("omegas", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omegaw_
    (
        IOobject
        (
            IOobject::groupName("omegaw", U.group()),
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
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool twophasekOmegaVeg<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        Clambda_.readIfPresent(this->coeffDict());
        betaOmega_.readIfPresent(this->coeffDict());
        alphaOmegaOmega_.readIfPresent(this->coeffDict());
        alphaKOmega_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void twophasekOmegaVeg<BasicTurbulenceModel>::correct()
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
    volScalarField magS = 2*magSqr(Sij);
    volScalarField G
    (
        this->GName(),
        nut*2*magSqr(Sij)
    );

    // Update omega and G at the wall
    omegas_.boundaryFieldRef().updateCoeffs();

    volTensorField Omij(-skew(fvc::grad(U)));
    volVectorField Gradk(fvc::grad(ks_));
    volVectorField Gradomega(fvc::grad(omegas_));
    volScalarField alphadCheck_(Gradk & Gradomega);

    const volScalarField CDkOmega
    (
    sigmad_*pos(0.15-alpha)*pos(alphadCheck_)*(alphadCheck_)/omegas_
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
            mag((Omij & Omij & Sij)/(pow((Cmu_*omegas_), 3)))
        );
    }

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegasEqn
    (
        fvm::ddt(omegas_)
      + fvm::div(phi, omegas_)
      - fvm::Sp(fvc::div(phi), omegas_)
      - fvm::laplacian(DomegaEff(), omegas_)
      ==
      - fvm::SuSp (-alphaOmega_*G/ks_, omegas_)
      //- fvm::SuSp (-alphaOmega_*magS/omegas_, omegas_)
      - fvm::Sp(ESD_, omegas_)
      - fvm::Sp
        (
            betaOmega_*
            (
                (scalar(1.0)+scalar(85.0)*XsiOmega())
                /(scalar(1.0)+scalar(100.0)*XsiOmega())
            )*omegas_(),
            omegas_
        )
      + CDkOmega
      + ESD2()*fvm::Sp(C3om_*KE2_, omegas_)
      //+ fvm::Sp((C4om_*KE4_*ESD5_/omegas_), omegas_)
      + fvm::Sp((C4om_*KE4_*ESD5_*nut/ks_), omegas_)
    );
    if (writeTke_)
    {
        #include "writeTKEBudget_kOmegaVeg.H"
    }
    omegasEqn.ref().relax();
    fvOptions.constrain(omegasEqn.ref());
    omegasEqn.ref().boundaryManipulate(omegas_.boundaryFieldRef());
    solve(omegasEqn);
    fvOptions.correct(omegas_);
    bound(omegas_, this->omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> ksEqn
    (
        fvm::ddt(ks_)
      + fvm::div(phi, ks_)
      - fvm::Sp(fvc::div(phi), ks_)
      - fvm::laplacian(DkEff(), ks_)
      ==
      - fvm::SuSp(-G/ks_, ks_)
      //- fvm::SuSp(-magS/omegas_, ks_)
      + fvm::Sp(-Cmu_*omegas_, ks_)
      + fvm::Sp(ESD_, ks_)
      //+ fvm::Sp(KE4_*ESD4_/omegas_, ks_)
      + fvm::Sp(KE4_*ESD4_*nut/ks_, ks_)
      + ESD2()*fvm::Sp(KE2_, ks_)
      //Spectral shortcut
      - fvm::Sp(KE8_*ESD8_, ks_)
    );
    
    ksEqn.ref().relax();
    fvOptions.constrain(ksEqn.ref());
    solve(ksEqn);
    fvOptions.correct(ks_);
    bound(ks_, this->kMin_);

    // Wake dissipation
    omegaw_ = KE7_*ESD7_/Clambda_*sqrt(kw_);
    fvOptions.correct(omegaw_);
    bound(omegaw_, this->omegaMin_);

    // Wake turbulent kinetic energy equation
    tmp<fvScalarMatrix> kwEqn
    (
    fvm::ddt(kw_)
      + fvm::div(phi, kw_)
      - fvm::Sp(fvc::div(phi), kw_)
      - fvm::laplacian(DkEff(), kw_)
      ==
      //turbulent production by vegetation drag force
      - fvm::SuSp(-KE6_*ESD6_/kw_, kw_)
      //turbulent dissipation by vegetation drag force
      + fvm::Sp(-Clambda_*omegaw_, kw_)
      //Spectral shortcut
      + fvc::Sp(KE8_*ESD8_, ks_)
    );
    
    kwEqn.ref().relax();
    fvOptions.constrain(kwEqn.ref());
    solve(kwEqn);
    fvOptions.correct(kw_);
    bound(kw_, this->kMin_);

    k_ = ks_ + kw_;
    omega_ = max(omegas_, omegaw_);

    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);

    fvOptions.correct(k_);
    bound(k_, this->kMin_);


    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
