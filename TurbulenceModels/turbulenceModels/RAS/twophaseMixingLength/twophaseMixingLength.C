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

\*---------------------------------------------------------------------------*/

#include "twophaseMixingLength.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void twophaseMixingLength<BasicTurbulenceModel>::correctNut()
{
    scalar cmu34 = pow(Cmu_.value(), 3.0/4.0);
    this->nut_ = 0.8*cmu34*k_*k_/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
twophaseMixingLength<BasicTurbulenceModel>::twophaseMixingLength
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
    eddyViscosity<RASModel<BasicTurbulenceModel>>
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
    twophaseRASProperties_
    (
        IOobject
        (
            "twophaseRASProperties",
            this->runTime_.constant(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    expoLM_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "expoLM",
            this->coeffDict_,
            1.0
        )
    ),
    alphaMaxLM_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaMaxLM",
            this->coeffDict_,
            0.55
        )
    ),
    kappaLM_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappaLM",
            this->coeffDict_,
            0.225
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
bool twophaseMixingLength<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
      return true;
    }
    else
    {
      return false;
    }
}
template<class BasicTurbulenceModel>
void twophaseMixingLength<BasicTurbulenceModel>::correct()
{
    if (not this->turbulence_)
    {
        return;
    }

//
// Mixing Leng turbulence model
// (only for 1D cases with Y the wall-normal direction)
//
    dimensionedScalar LmSmall
    (
        "LmSmall",
        dimensionSet(0, 1, 0, 0, 0, 0, 0),
        scalar(1e-4)
    );

// Local references
    const volScalarField& alpha = this->alpha_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> gradU = fvc::grad(U);
    volSymmTensorField D = symm(gradU);
    volScalarField magD = ::sqrt(2.0)*mag(D);
    gradU.clear();

    volVectorField centres = U.mesh().C();
    volScalarField Y = centres.component(1);

    scalar Lm = 0.;
    scalar dY;
    scalar expoLM=expoLM_.value();
    scalar kappaLMs=kappaLM_.value();
    scalar alphaMaxLMs = alphaMaxLM_.value();
    scalar LmPhi = 0.;
    scalar cmu34 = pow(Cmu_.value(), 3.0/4.0);


    nut.storePrevIter();
    forAll(U, cellI)
    {
        if (cellI==0)
        {
            dY = Y[cellI+1]-Y[cellI];
        }
        else
        {
            dY = Y[cellI]-Y[cellI-1];
        }
        LmPhi = LmPhi
              + kappaLMs*max
              (
                  scalar(1.0)
                - Foam::pow(min(alpha[cellI]/alphaMaxLMs, scalar(1.0)),
                            expoLM), 0.
              )*dY;
        Lm = LmPhi;
        nut[cellI] = pow(Lm, 2)*magD[cellI];
        // for kinetic theory k is required
        k_[cellI] = pow(nut[cellI]/(0.8*max(Lm, 1e-4)), 2);
        epsilon_[cellI] = cmu34*Foam::pow(k_[cellI], 1.5)/max(Lm, 1e-4);
    }
    nut.relax();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
