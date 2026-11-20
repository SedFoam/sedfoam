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

#include "twophasekEpsilonVeg.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void twophasekEpsilonVeg<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(ks_)/epsilons_ + Clambda_*sqr(kw_)/epsilonw_;
    this->nut_.min(nutMax_);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
twophasekEpsilonVeg<BasicTurbulenceModel>::twophasekEpsilonVeg
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
    writeTke_
    (
        Switch::getOrAddToDict
        (
            "writeTke",
            this->coeffDict_,
            false
        )
    ),
    C3ep_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C3ep",
            this->coeffDict_,
            1.2
        )
    ),
    C4ep_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C4ep",
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
    C1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
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
    alphak_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphak",
            this->coeffDict_,
            1.0
        )
    ),
    alphaEps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaEps",
            this->coeffDict_,
            1.3
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
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilons_
    (
        IOobject
        (
            IOobject::groupName("epsilons", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilonw_
    (
        IOobject
        (
            IOobject::groupName("epsilonw", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool twophasekEpsilonVeg<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict_);
        Clambda_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict_);
        C2_.readIfPresent(this->coeffDict_);
        alphak_.readIfPresent(this->coeffDict_);
        alphaEps_.readIfPresent(this->coeffDict_);
        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void twophasekEpsilonVeg<BasicTurbulenceModel>::correct()
{
    if (not this->turbulence_)
    {
        return;
    }

    // Local references
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const surfaceScalarField& phi = this->phi_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    volTensorField GradU(fvc::grad(U));
    volSymmTensorField Sij(symm(GradU));

    //volScalarField::Internal G
    volScalarField G
    (
        this->GName(),
        nut*(GradU && dev(twoSymm(GradU)))
    );

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsilonsEqn
    (
        fvm::ddt(epsilons_)
      + fvm::div(phi, epsilons_)
      - fvm::Sp(fvc::div(phi), epsilons_)
      - fvm::laplacian(DepsilonEff(), epsilons_)
     ==
      - fvm::SuSp(-C1_*G/ks_, epsilons_)
      + fvm::Sp(-C2_*epsilons_/ks_, epsilons_)
      + fvm::Sp(C3ep_*ESD_, epsilons_)
      + ESD2()*fvm::Sp(C3ep_*KE2_, epsilons_)
      + fvm::Sp(C4ep_*KE4_*ESD5_*nut/ks_, epsilons_)
    );
    if (writeTke_)
    {
        #include "writeTKEBudget_kEpsilonVeg.H"
    }
    epsilonsEqn.ref().relax();
    fvOptions.constrain(epsilonsEqn.ref());
    epsilonsEqn.ref().boundaryManipulate(epsilons_.boundaryFieldRef());
    solve(epsilonsEqn);
    fvOptions.correct(epsilons_);
    bound(epsilons_, this->epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> ksEqn
    (
        fvm::ddt(ks_)
      + fvm::div(phi, ks_)
      - fvm::Sp(fvc::div(phi), ks_)
      - fvm::laplacian(DkEff(), ks_)
     ==
      - fvm::SuSp(-G/ks_, ks_)
      - fvm::Sp(epsilons_/ks_, ks_)
      + fvm::Sp(ESD_, ks_)
      + fvm::Sp(KE4_*ESD4_*nut/ks_, ks_)
      + ESD2()*fvm::Sp(KE2_, ks_)
      //Spectral shortcut
      - fvc::Sp(KE8_*ESD8_, ks_)
    );

    ksEqn.ref().relax();
    fvOptions.constrain(ksEqn.ref());
    solve(ksEqn);
    fvOptions.correct(ks_);
    bound(ks_, this->kMin_);

    // Wake dissipation
    epsilonw_ = KE7_*ESD7_*kw_*sqrt(kw_);
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
      + fvm::Sp(-epsilonw_/kw_, kw_)
      //Spectral shortcut
      + fvc::Sp(KE8_*ESD8_, ks_)
    );

    kwEqn.ref().relax();
    fvOptions.constrain(kwEqn.ref());
    solve(kwEqn);
    fvOptions.correct(kw_);
    bound(kw_, this->kMin_);

    k_ = ks_ + kw_;
    epsilon_ = epsilons_ + epsilonw_;

    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
