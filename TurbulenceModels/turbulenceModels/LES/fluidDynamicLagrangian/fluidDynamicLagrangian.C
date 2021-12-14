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

#include "fluidDynamicLagrangian.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void fluidDynamicLagrangian<BasicTurbulenceModel>::correctNut
(
    const tmp<volTensorField>& gradU
)
{
    this->nut_ = min((flmb_/fmmb_)*sqr(this->delta())*mag(dev(symm(gradU))),
                     nutMax_);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void fluidDynamicLagrangian<BasicTurbulenceModel>::correctNut()
{
    correctNut(fvc::grad(this->U_));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
fluidDynamicLagrangian<BasicTurbulenceModel>::fluidDynamicLagrangian
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    flmb_
    (
        IOobject
        (
            IOobject::groupName("flmb", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    fmmb_
    (
        IOobject
        (
            IOobject::groupName("fmmb", this->U_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    Cyb_
    (
        IOobject
        (
            "Cyb",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),
    theta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "theta",
            this->coeffDict_,
            1.5
        )
    ),
    D_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "D",
            this->coeffDict_,
            0.0
        )
    ),
    nutMax_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "nutbMax",
            this->coeffDict_,
            dimLength*dimVelocity,
            1e-3
        )
    ),
    CybMax_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CybMax",
            this->coeffDict_,
            1.0
        )
    ),
    simpleFilter_(U.mesh()),
    filterPtr_(LESfilter::New(U.mesh(), this->coeffDict())),
    filter_(filterPtr_()),

    flmb0_("flmb0", flmb_.dimensions(), 0.0),
    fmmb0_("fmmb0", fmmb_.dimensions(), VSMALL)
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool fluidDynamicLagrangian<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        filter_.read(this->coeffDict());
        theta_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void fluidDynamicLagrangian<BasicTurbulenceModel>::correct()
{
    if (not this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField S(dev(symm(gradU)));
    volScalarField magS(mag(S));
    volScalarField magSqrS(magSqr(S));

    volVectorField Uf(filter_(U));
    volSymmTensorField Sf(dev(symm(fvc::grad(Uf))));
    volScalarField magSf(mag(Sf));
    volScalarField magSqrSf(magSqr(Sf));

    volSymmTensorField L(dev(filter_(sqr(U)) - (sqr(filter_(U)))));
    volSymmTensorField M
    (
        2.0*sqr(this->delta())*(filter_(magS*S) - 4.0*magSf*Sf)
    );

    volScalarField invT
    (
        alpha*rho*(1.0/(theta_.value()*this->delta()))
        *pow(flmb_*fmmb_, 1.0/8.0)
    );

    volScalarField LM(L && M);

    fvScalarMatrix flmEqn
    (
        fvm::ddt(alpha, rho, flmb_)
      + fvm::div(alphaRhoPhi, flmb_)
     ==
        invT*LM
      - fvm::Sp(invT, flmb_)
      + fvOptions(alpha, rho, flmb_)
    );

    flmEqn.relax();
    fvOptions.constrain(flmEqn);
    flmEqn.solve();
    fvOptions.correct(flmb_);
    bound(flmb_, flmb0_);

    volScalarField MM(M && M);

    fvScalarMatrix fmmEqn
    (
        fvm::ddt(alpha, rho, fmmb_)
      + fvm::div(alphaRhoPhi, fmmb_)
     ==
        invT*MM
      - fvm::Sp(invT, fmmb_)
      + fvOptions(alpha, rho, fmmb_)
    );

    fmmEqn.relax();
    fvOptions.constrain(fmmEqn);
    fmmEqn.solve();
    fvOptions.correct(fmmb_);
    bound(fmmb_, fmmb0_);

    correctNut(gradU);
    
    volScalarField Lkk(1.0/3.0*tr((filter_(sqr(U)))-sqr(filter_(U))));

    volScalarField Mkk
    (
        -2.0*sqr(this->delta())*(filter_(magSqrS) - 4.0*magSqrSf)
    );

    Lkk.correctBoundaryConditions();
    Mkk.correctBoundaryConditions();

    volScalarField Lkka(fvc::average(Lkk*Mkk));
    volScalarField Mkka(fvc::average(Mkk*Mkk));

    Cyb_ = min(Lkka/Mkka, CybMax_);
    bound(Cyb_, VSMALL);

}

template<class BasicTurbulenceModel>
tmp<volScalarField> fluidDynamicLagrangian<BasicTurbulenceModel>::
spherSigmaSGS()
{
    const volVectorField& U = this->U_;
    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField S(dev(symm(gradU)));
    volScalarField magSqrS(magSqr(S));

    const tmp<volScalarField> spherSigma(Cyb_*sqr(this->delta())*magSqrS);

    return spherSigma;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> fluidDynamicLagrangian<BasicTurbulenceModel>::
sqrDeltaD()
{
    const tmp<volScalarField> sqrDeltaD(sqr(this->delta())*D_);
    return sqrDeltaD;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
