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
    const Foam::volScalarField& pa,
    const Foam::dimensionedScalar& Dsmall
)
:
    alpha_(phasea.alpha()),
    phia_(phasea.phi()),
    rhoa_(phasea.rho()),
    da_(phasea.d()),
    rhob_(phaseb.rho()),
    nub_(phaseb.nu()),

    pa_new_value(pa),

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
    granularRheology_
    (
        granularRheologyProperties_.get<Switch>("granularRheology")
    ),
    granularDilatancy_
    (
        granularRheologyProperties_.get<Switch>("granularDilatancy")
    ),
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
    alphaMaxG_
    (
        granularRheologyProperties_.getOrDefault
        (
            "alphaMaxG",
            dimensionedScalar("alphaMaxG",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.6)
        )
    ),
    mus_
    (
        granularRheologyProperties_.getOrDefault
        (
            "mus",
            dimensionedScalar("mus",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.38)
        )
    ),
    mu2_
    (
        granularRheologyProperties_.getOrDefault
        (
            "mu2",
            dimensionedScalar("mu2",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.64)
        )
    ),
    I0_
    (
        granularRheologyProperties_.getOrDefault
        (
            "I0",
            dimensionedScalar("I0",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.3)
        )
    ),
    Bphi_
    (
        granularRheologyProperties_.getOrDefault
        (
            "Bphi",
            dimensionedScalar("Bphi",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.31)
        )
    ),
    n_
    (
        granularRheologyProperties_.getOrDefault
        (
            "n",
            dimensionedScalar("n",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          2.5)
        )
    ),
    BulkFactor_
    (
        granularRheologyProperties_.getOrDefault
        (
                "BulkFactor",
                dimensionedScalar("BulkFactor",
                    dimensionSet(0, 0, 0, 0, 0, 0, 0),
                    0)
        )
    ),
    alpha_c_
    (
        granularRheologyProperties_.getOrDefault
        (
                "alpha_c",
                dimensionedScalar("alpha_c",
                    dimensionSet(0, 0, 0, 0, 0, 0, 0),
                    0.585)
        )
    ),

    K_dila_
    (
        granularRheologyProperties_.getOrDefault
        (
                "K_dila",
                dimensionedScalar("K_dila",
                    dimensionSet(0, 0, 0, 0, 0, 0, 0),
                    0)
        )
    ),
    relaxPa_
    (
        granularRheologyProperties_.getOrDefault
        (
            "relaxPa",
            dimensionedScalar("relaxPa",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          1)
        )
    ),
    PaMin_
    (
        granularRheologyProperties_.getOrDefault
        (
            "PaMin",
            dimensionedScalar("PaMin",
                          dimensionSet(0, 1, 0, 0, 0, 0, 0),
                          1e-6)
        )
    ),
    tau_inv_min_
    (
        granularRheologyProperties_.getOrDefault
        (
            "tau_inv_min",
            dimensionedScalar("tau_inv_min",
                          dimensionSet(0, 0, -1, 0, 0, 0, 0),
                          1e-12)
        )
    ),
    muI_
    (
        IOobject
        (
            "muI",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("zero", alpha_.dimensions(), 0.0)
    ),
    mua_
    (
        IOobject
        (
            "mua",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    p_p_total_
    (
        IOobject
        (
            "p_p_total",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),

    delta_
    (
        IOobject
        (
            "delta",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("zero", alpha_.dimensions(), 0.0)
    ),

    nuvb_
    (
        IOobject
        (
            "nuvb",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),
    I_
    (
        IOobject
        (
            "I",
            alpha_.time().timeName(),
            alpha_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::granularRheologyModel::~granularRheologyModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::granularRheologyModel::solve
(
    const volScalarField& magD,
    const volScalarField& pf,
    const dimensionedScalar& alphaSmall,
    const dimensionedScalar& Dsmall
)
{
    if (granularRheology_ == false)
    {
        return;
    }
    dimensionedScalar Dsmall2
    (
        "Dsmall2",
        dimensionSet(0, 0, -2, 0, 0, 0, 0),
        1e-8
    );
    Dsmall2 = sqr(Dsmall);
    //
    // compute the particulate velocity shear rate
    //
    volScalarField magD2(pow(magD, 2));

    //
    // Shear induced particulate pressure
    //
    pa_ = PPressureModel_->pa
    (
        pf, Bphi_, rhoa_, da_, rhob_, nub_, magD,
        alpha_, alphaMaxG_, alphaSmall
    );

    // Relaxing shear induced particulate pressure
    //  relaxPa_ controls the relaxation of pa. Low values lead to relaxed pa
    //  whereas large value are prone to numerical error
    volScalarField tau_inv_par(relaxPa_*alpha_*magD);
    tau_inv_par.max(tau_inv_min_);

    fvScalarMatrix paEqn
    (
         fvm::ddt(pa_new_value)
         + fvm::div(phia_, pa_new_value, "div(phia,pa_new_value)")
         - fvm::Sp(fvc::div(phia_), pa_new_value)
        ==
        tau_inv_par*(pa_)
        -fvm::Sp(tau_inv_par, pa_new_value)
    );
    paEqn.relax();
    paEqn.solve();

    pa_new_value.max(0.0);

    pa_=pa_new_value;
    pa_.correctBoundaryConditions();
    //total particle pressure(shear induced+contact contributions)
    p_p_total_ = mag(pa_new_value+pf);
    p_p_total_.max(PaMin_);

    //  Compute the particulate friction coefficient
    muI_ = FrictionModel_->muI(mus_, mu2_, I0_, p_p_total_, rhoa_, da_, rhob_,
                               nub_, magD, Dsmall);

    //  Compute the inertial/viscous number
    I_ = FrictionModel_->I(p_p_total_, rhoa_, da_, rhob_, nub_, magD);

    // Dilatancy model
    if (granularDilatancy_)
    {
        volScalarField alphaEq_
        (
            PPressureModel_->alphaEq
            (
                p_p_total_,
                Bphi_,
                rhoa_,
                da_,
                rhob_,
                nub_,
                magD,
                alpha_c_
            )
        );
        delta_ = K_dila_*(alpha_ - alphaEq_);

        delta_.min( 0.4);
        delta_.max(-0.4);
    }
    //  Compute the regularized particulate viscosity
    mua_ = muI_* p_p_total_ / pow(magD2 + Dsmall2, 0.5);

    // Compute bulk viscosity (by default BulkFactor = 0)s
    lambda_ = BulkFactor_*p_p_total_ / pow(magD2 + Dsmall2, 0.5);

    // Compute the Effective fluid viscosity
    nuvb_ = FluidViscosityModel_->nuvb(alpha_, nub_, alphaMaxG_, alphaSmall,
                                       n_);
}
// ************************************************************************* //
