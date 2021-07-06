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

#include "twophaseMixingLengthConst.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void twophaseMixingLengthConst<BasicTurbulenceModel>::correctNut()
{
    scalar cmu34 = 0.1643;
    this->nut_ = 0.8*cmu34*k_*k_/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
twophaseMixingLengthConst<BasicTurbulenceModel>::twophaseMixingLengthConst
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
bool twophaseMixingLengthConst<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void twophaseMixingLengthConst<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
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
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    volScalarField magD(::sqrt(2.0)*mag(symm(fvc::grad(U))));

    volScalarField Y(U.mesh().C().component(vector::Y));

    scalar Lm(0.);

    scalar cmu34 = 0.1643;


    nut.storePrevIter();
    forAll(U, cellI)
    {
        Lm = 0.2*0.007;
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
