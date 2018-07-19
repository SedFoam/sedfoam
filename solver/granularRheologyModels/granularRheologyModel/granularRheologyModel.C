/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "granularRheologyModel.H"
#include "FrictionModel.H"
#include "PPressureModel.H"
#include "FluidViscosityModel.H"
#include "surfaceInterpolate.H"
#include "mathematicalConstants.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::granularRheologyModel::granularRheologyModel
(
    const Foam::phaseModel& phasea,
    const Foam::phaseModel& phaseb,
    const Foam::volScalarField& pa
)
:
    alpha_(phasea.alpha()),
    phia_(phasea.phi()),
    paOld(pa),

    rhoa_(phasea.rho()),
    da_(phasea.d()),

    rhob_(phaseb.rho()),
    nub_(phaseb.nu()),

    granularRheologyProperties_
    (
        IOobject
        (
            "granularRheologyProperties",
            alpha_.time().constant(),
            alpha_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    granularRheology_(granularRheologyProperties_.lookup("granularRheology")),
    FrictionModel_
    (
        granularRheologyModels::FrictionModel::New
        (
            granularRheologyProperties_
        )
    ),
    PPressureModel_
    (
        granularRheologyModels::PPressureModel::New
        (
            granularRheologyProperties_
        )
    ),
    FluidViscosityModel_
    (
        granularRheologyModels::FluidViscosityModel::New
        (
            granularRheologyProperties_
        )
    ),
    alphaMaxG_(granularRheologyProperties_.lookup("alphaMaxG")),
    mus_(granularRheologyProperties_.lookup("mus")),
    mu2_(granularRheologyProperties_.lookup("mu2")),
    I0_(granularRheologyProperties_.lookup("I0")),
    Bphi_(granularRheologyProperties_.lookup("Bphi")),
    n_(granularRheologyProperties_.lookup("n")),
    Dsmall_(granularRheologyProperties_.lookup("Dsmall")),
    relaxPa_(granularRheologyProperties_.lookup("relaxPa")),
    muI_
    (
        IOobject
        (
            "muI",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),
    mua_
    (
        IOobject
        (
            "mua",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),

    pa_
    (
        IOobject
        (
            "pa",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),

     nuvb_
     (
        IOobject
        (
            "nuvb",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
     )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularRheologyModel::~granularRheologyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::granularRheologyModel::solve
(
    const volTensorField& gradUat,
    const volScalarField& pf,
    const dimensionedScalar& alphaSmall,
    const dimensionedScalar& dt
)
{
    if (!granularRheology_)
    {
        return;
    }
    dimensionedScalar Dsmall2
    (
        "Dsmall2",
        dimensionSet(0 , 0 ,-2 ,0 , 0, 0, 0),
        1e-8
    );
    Dsmall2 = sqr(Dsmall_);
    //
    // compute the particulate velocity shear rate
    //
    //volSymmTensorField D = dev(symm(gradUat));
    volSymmTensorField D = symm(gradUat);
   
    volScalarField magD = ::sqrt(2.0)*mag(D);
    volScalarField magD2 = pow(magD,2);

    volScalarField patot_ = pf*scalar(0.0);

    //
    //      Shear induced particulate pressure 
    //
    pa_ = PPressureModel_->pa
    (
        pf, Bphi_, rhoa_, da_, rhob_, nub_, magD,
        alpha_, alphaMaxG_, alphaSmall
    );
    // relaxation of the shear induced particle pressure
    pa_ = paOld + relaxPa_*(pa_-paOld);

    // Add contact pressure to the shear induced contribution
    patot_ = pa_ + pf;    
    
    //  Compute the particulate friction coefficient
    muI_ = FrictionModel_->muI(mus_, mu2_, I0_, patot_, rhoa_, da_, rhob_, nub_, magD, Dsmall_);

    //  Compute the regularized particulate viscosity
    mua_ = muI_* patot_ / pow(magD2 + Dsmall2, 0.5);
    //mua_ = muI_ * patot_ / (magD+Dsmall_);
   
// Compute mua_ on the bottom patch
    forAll(alpha_.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(alpha_.boundaryField()[patchi]))
        {
            mua_.boundaryFieldRef()[patchi] =
            (
                muI_.boundaryFieldRef()[patchi]*patot_.boundaryFieldRef()[patchi]
                /pow(magD2.boundaryFieldRef()[patchi] + Dsmall2.value(), 0.5)
            );
        }
/*        if (isA<zeroGradientFvPatchVectorField>(Ua_.boundaryField()[patchi]))
        {
            mua_.boundaryFieldRef()[patchi] = 0;
            muI_.boundaryFieldRef()[patchi] = 0;
            //Info <<mua_.boundaryFieldRef()[patchi]<<nl<<endl;
        }
*/
    }
    // Set bulk viscosity to zero
    lambda_ = scalar(0.0)*mua_;

    // Compute the Effective fluid viscosity
    nuvb_ = FluidViscosityModel_->nuvb(alpha_, nub_, alphaMaxG_, alphaSmall, n_);
}
// ************************************************************************* //
