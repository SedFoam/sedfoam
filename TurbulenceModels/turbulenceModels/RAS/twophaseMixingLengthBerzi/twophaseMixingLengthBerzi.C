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

#include "twophaseMixingLengthBerzi.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void twophaseMixingLengthBerzi<BasicTurbulenceModel>::correctNut()
{
    scalar cmu34 = pow(Cmu_.value(), 3.0/4.0);
    this->nut_ = 0.8*cmu34*k_*k_/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
twophaseMixingLengthBerzi<BasicTurbulenceModel>::twophaseMixingLengthBerzi
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
    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    expoLM_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "expoLM",
            this->coeffDict_,
            1.0
        )
    ),
    alphaMaxLM_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaMaxLM",
            this->coeffDict_,
            0.55
        )
    ),
    kappaLM_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappaLM",
            this->coeffDict_,
            0.225
        )
    ),
    d_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "d",
            this->coeffDict_,
            0.005
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
    )
{

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool twophaseMixingLengthBerzi<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        return true;
    }

    return false;
}

template<class BasicTurbulenceModel>
void twophaseMixingLengthBerzi<BasicTurbulenceModel>::correct()
{
    if (not this->turbulence_)
    {
        return;
    }

//
// Mixing Length turbulence model
// (only for 1D cases with Y the wall-normal direction)
//
    dimensionedScalar LmSmall
    (
        "LmSmall",
        dimensionSet(0, 1, 0, 0, 0, 0, 0),
        scalar(1e-4)
    );

// Local references
    const volScalarField& beta = this->alpha_;
    const volScalarField alpha = 1 - beta;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    volScalarField magD(::sqrt(2.0)*mag(symm(fvc::grad(U))));

    volScalarField Y(U.mesh().C().component(vector::Y));

    scalar Lm(0.);
    scalar dY;
    scalar expoLM(expoLM_.value());
    scalar kappaLMs(kappaLM_.value());
    scalar d(d_.value());
    scalar alphaMaxLMs(alphaMaxLM_.value());
    scalar LmPhi(0.);
    scalar cmu34(pow(Cmu_.value(), 3.0/4.0));

    nut.storePrevIter();
    forAll(U, cellI)
    {
	Lm = max(3*d*pow(alphaMaxLMs-alpha[cellI],expoLM), 0.2*d);
        nut[cellI] = pow(Lm, 2)*magD[cellI];
        // for kinetic theory k is required
        k_[cellI] = pow(nut[cellI]/(0.1*max(Lm, 1e-4)), 2);
        epsilon_[cellI] = cmu34*Foam::pow(k_[cellI], 1.5)/max(Lm, 1e-4);
    }
    nut.relax();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
