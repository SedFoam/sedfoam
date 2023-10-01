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
#include "kineticTheoryModel.H"
#include "surfaceInterpolate.H"
#include "mathematicalConstants.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::kineticTheoryModel
(
    const Foam::phaseModel& phasea,
    const Foam::volVectorField& Ub,
    const Foam::dragModel& draga
)
:
    phasea_(phasea),
    Ua_(phasea.U()),
    Ub_(Ub),
    alpha_(phasea.alpha()),
    phia_(phasea.phi()),
    draga_(draga),
    rhoa_(phasea.rho()),
    da_(phasea.d()),
    nua_(phasea.nu()),

    kineticTheoryProperties_
    (
        IOobject
        (
            "kineticTheoryProperties",
            Ua_.time().constant(),
            Ua_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    kineticTheory_
    (
        kineticTheoryProperties_.get<Switch>("kineticTheory")
    ),
    writeTBudget_
    (
        kineticTheoryProperties_.getOrDefault<Switch>("writeTBudget", false)
    ),
    extended_
    (
        kineticTheoryProperties_.getOrDefault<Switch>("extended", false)
    ),
    limitProduction_
    (
        kineticTheoryProperties_.getOrDefault<Switch>("limitProduction", false)
    ),
    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    conductivityModel_
    (
        conductivityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    pseudoConductivityModel_
    (
        pseudoConductivityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    radialModel_
    (
        radialModel::New
        (
            kineticTheoryProperties_
        )
    ),
    granularPressureModel_
    (
        granularPressureModel::New
        (
            kineticTheoryProperties_
        )
    ),
    saltationModel_
    (
        saltationModel::New
        (
            kineticTheoryProperties_
        )
    ),
    e_
    (
        kineticTheoryProperties_.getOrDefault
        (
            "e",
            dimensionedScalar("e",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.9)
        )
    ),
    alphaMax_
    (
        kineticTheoryProperties_.getOrDefault
        (
            "alphaMax",
            dimensionedScalar("alphaMax",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.6)
        )
    ),
    MaxTheta
    (
        kineticTheoryProperties_.getOrDefault
        (
            "MaxTheta",
            dimensionedScalar("MaxTheta",
                          dimensionSet(0, 2, -2, 0, 0, 0, 0),
                          1e3)
        )
    ),
    phi_
    (
        kineticTheoryProperties_.getOrDefault
        (
            "phi",
            dimensionedScalar("phi",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          32)
        )*M_PI/180.0 //32° angle of repose
    ),
    muPart_
    (
        kineticTheoryProperties_.getOrDefault
        (
            "muPart",
            dimensionedScalar("muPart",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.0)
        )
    ),
    killJ2_
    (
        kineticTheoryProperties_.getOrDefault
        (
            "killJ2_",
            dimensionedScalar
            (
                "killJ2_",
                dimensionSet(0, 0, 0, 0, 0, 0, 0),
                1
            )
        )
    ),
    quadraticCorrectionJ1_
    (
        kineticTheoryProperties_.getOrDefault
        (
            "quadraticCorrectionJ1_",
            dimensionedScalar
            (
                "quadraticCorrectionJ1_",
                dimensionSet(0, 0, 0, 0, 0, 0, 0),
                1
            )
        )
    ),
    Theta_
    (
        IOobject
        (
            "Theta",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0), 0.0)
    ),
    mua_
    (
        IOobject
        (
            "mua",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    muf_
    (
        IOobject
        (
            "muf",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    muSaltCoef_
    (
        IOobject
        (
            "muSaltCoef",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    pa_
    (
        IOobject
        (
            "pa",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    pf_
    (
        IOobject
        (
            "pf_",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    ppMagf_
    (
        IOobject
        (
            "ppMagf",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    kappaSaltCoef_
    (
        IOobject
        (
            "kappaSaltCoef",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),
    kappaAlpha_
    (
        IOobject
        (
            "kappaAlpha",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, 1, -3, 0, 0), 0.0)
    ),
    gs0_
    (
        IOobject
        (
            "gs0",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    gs0Prime_
    (
        IOobject
        (
            "gs0prime",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::~kineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kineticTheoryModel::solve
(
    const volTensorField& gradUat,
    const volScalarField& kb,
    const volScalarField& epsilonb,
    const volScalarField& nutf,
    const dimensionedScalar& B,
    const dimensionedScalar& tt
)
{
    if (kineticTheory_ == false)
    {
        return;
    }
    ////////////////////////////////
    //Define some usefull quantities
    ////////////////////////////////
    const scalar sqrtPi(sqrt(constant::mathematical::pi));
    const scalar Pi(constant::mathematical::pi);
    dimensionedScalar alphaSmall
    (
        "small",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        scalar(1.0e-4)
    );
    dimensionedScalar Tsmall
    (
        "small",
        dimensionSet(0, 2, -2, 0, 0, 0, 0),
        SMALL
    );
    dimensionedScalar Ksmall
    (
        "small",
        dimensionSet(1, -3, -1, 0, 0, 0, 0),
        scalar(1.0e-4)
    );
    volScalarField ThetaSqrt(sqrt(Theta_));
    dimensionedScalar Tpsmall_
    (
         "Tpsmall_",
         dimensionSet(1, -1, -3, 0, 0, 0, 0),
         scalar(1e-10)
    );
    volTensorField dU(fvc::grad(Ua_));//gradUat.T(); //that is fvc::grad(Ua_);
    volSymmTensorField D(symm(dU));   //0.5*(dU + dU.T)

    volScalarField ThetaClip(Theta_);
    if (limitProduction_) // limit the production by clipping Theta to 2/3kb
    {
        ThetaClip = min(Theta_, (2.0/3.0)*max(kb));
    }

    //////////////////////////////////
    // Radial distribution function g0
    //////////////////////////////////
    gs0_ = radialModel_->g0(min(alpha_, alphaMax_-alphaSmall),
                            alphaMax_, muPart_);
    gs0Prime_ = radialModel_->g0prime(min(alpha_, alphaMax_ - alphaSmall),
                                      alphaMax_, muPart_);

    //////////////////////////////////////////
    // collisional pressure  (Eq. 3.22, p. 45)
    // ///////////////////////////////////////
    volScalarField PsCoeff
    (
        granularPressureModel_->granularPressureCoeff
        (
            alpha_,
            gs0_,
            rhoa_,
            e_
        )
    );

    ///////////////////////////////////////
    // Saltation Model
    ///////////////////////////////////////
    volScalarField K(draga_.K(mag(Ua_ - Ub_)));
    muSaltCoef_ = saltationModel_->musalt(alpha_+alphaSmall, ThetaClip+Tsmall,
                    rhoa_, da_, K+Ksmall);
    kappaSaltCoef_ = saltationModel_->kappasalt(alpha_+alphaSmall,
                    ThetaClip+Tsmall, rhoa_, da_, K+Ksmall);

    ///////////////////////////////////////
    // Granular viscosity (Table 3.2, p.47)
    // and shear stress
    ///////////////////////////////////////
    mua_ = viscosityModel_->mua(alpha_+alphaSmall, ThetaClip+Tsmall,
                    gs0_, muSaltCoef_, K, rhoa_, da_, e_);
    // bulk viscosity  p. 45.
    lambda_ = viscosityModel_->lambda(alpha_+alphaSmall, ThetaClip+Tsmall,
                    gs0_, rhoa_, da_, e_);

    // stress tensor, Definitions, Table 3.1, p. 43
    volSymmTensorField tau
    (
        2.0*mua_*D + (lambda_ - (2.0/3.0)*mua_)*tr(D)*I
    );

    //////////////////////////////////////////////
    // Temperature conductivity (Table 3.3, p. 49)
    //////////////////////////////////////////////
    kappa_ = conductivityModel_->kappa(alpha_, ThetaClip, gs0_, kappaSaltCoef_,
                    K, rhoa_, da_, e_);

    //////////////////////////////////////////////
    // Temperature conductivity (Table 3.3, p. 49)
    //////////////////////////////////////////////
    kappaAlpha_ = pseudoConductivityModel_->kappaAlpha(alpha_, ThetaClip, gs0_,
                    gs0Prime_, rhoa_, da_, e_);

    //////////////////////
    // Contact dissipation
    //////////////////////

    // Correlation length for the extended kinetic theory
    volScalarField Lc
    (
        da_*(alpha_+alphaSmall)/(alpha_+alphaSmall)
    );
    if (extended_)
    {// extended kinetic theory Jenkins (2007)
        Lc = da_*max
        (
            scalar(1),
            //0.5*pow(30./(1.+sqr(sqrtPi)/12.)*(1-e_)*sqr(alpha_)*gs0_, 1./3.)
            //0.5*pow(30*sqr(alpha_+alphaSmall)*(1-e_)/(1+sqr(sqrtPi)*(1+e_)*
            // (1-3*e_)/(1-1./4.*(1-sqr(e_))-5./24.*sqr(1-e_))), 1./3.)
            0.5*pow(30*pow(alpha_+alphaSmall, 2)*(1-e_)/
                     (1-Pi*(1+e_)*(1-3*e_)/
                     (1-1./4.*pow(1-e_, 2)-5./24.*(1-pow(e_, 2)))
                     ), 1./3.)
            *pow((alpha_+alphaSmall)*gs0_, 2./9.)
        );
    }
    // Inelastic dissipation (Eq. 3.24, p.50)
    volScalarField gammaCoeff
    (
        (3.0*(1.0 - sqr(e_))*sqr(alpha_)*rhoa_*gs0_
            *(4.0/Lc*ThetaSqrt/sqrtPi-tr(D)))
    );
    // Frictional dissipation (Chialvo and Sundaresan (2013) eqs. 22-23)
    //dimensionedScalar f_mu(3./2*muPart_*exp(-3*muPart_));
    // Frictional dissipation (Jenkins & Zhang (2002))
    dimensionedScalar f_mu(0);
    if (muPart_>alphaSmall)
    {
        //Tangential restitution coefficient assumed to be 0
        dimensionedScalar mu0 = 7./2*muPart_*(1+e_) + alphaSmall;
        dimensionedScalar a1 = muPart_/mu0*(Pi*mu0*(1-2/Pi*atan(mu0)) +
                     2*pow(mu0, 2)/(1+pow(mu0, 2))*(1-2*muPart_/mu0));
        dimensionedScalar a2 = 5./2*muPart_/mu0*(Pi/2*mu0*(1-2/Pi*atan(mu0)) +
                     (pow(mu0, 2)-pow(mu0, 4))/pow(1+pow(mu0, 2), 2));
        dimensionedScalar b1 = pow(muPart_, 2)/(1+pow(mu0, 2));
        dimensionedScalar b2 = 1./2*muPart_/mu0*(Pi/2*mu0*(1-2/Pi*atan(mu0)) +
                     pow(mu0, 2)/(1+pow(mu0, 2)));
        f_mu = 1./2*(a1-a2*b1/b2);
    }
    dimensionedScalar e_eff(e_ - f_mu);
    dimensionedScalar fric_correction((1-pow(e_eff, 2))/(1-pow(e_, 2)));
    gammaCoeff *= fric_correction;

    //////////////////////////////////////////////////////////
    // Work of the drag force by granular fluctuating velocity
    // NB, drag = K*alpha*beta,
    // this is inconsistent with momentum equation,
    // since the following form missed the drift velocity
    /////////////////////////////////////////////////////////
    // The following is to calculate parameter tmf_ in u_f*u_s correlation
    // (Danon et al. (1977))
    volScalarField flucVelCor_
    (
        Foam::exp(-B*rhoa_*6*epsilonb/max(kb*(1-alpha_)*K, Tpsmall_))
    );
    // Eq. 3.25, p. 50 J_int = J2 - J1*Theta_
    volScalarField J1(alpha_*(1-alpha_)*K*(3+2*quadraticCorrectionJ1_));
    volScalarField J2(alpha_*(1-alpha_)*K*(2*flucVelCor_*kb));

    ////////////////////////////////
    // Granular temperature equation
    // /////////////////////////////
    surfaceScalarField phi
    (
        1.5*rhoa_*phia_*fvc::interpolate((alpha_+alphaSmall))
    );
    // construct the granular temperature equation (Eq. 3.20, p. 44)
    // NB. note that there are two typos in Eq. 3.20
    // no grad in front of Ps
    // wrong sign in front of laplacian
    fvScalarMatrix ThetaEqn
    (
       fvm::ddt(1.5*(alpha_+alphaSmall)*rhoa_, Theta_)
     + fvm::div(phi, Theta_, "div(phi,Theta)")
     ==
       // Ps term.
       fvm::SuSp(-PsCoeff*fvc::div(phia_), Theta_)
       // production due to shear.
       + (tau && dU)
       // granular temperature conduction.
       + fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
       // granular temperature pseudo conduction.
       + fvc::laplacian(kappaAlpha_, alpha_, "laplacian(kappaAlpha,alpha)")
       // energy disipation due to inelastic collision.
       + fvm::Sp(-gammaCoeff, Theta_)
       // dissipation through drag force
       + fvm::Sp(-J1, Theta_)
       // transfer of fluctuating kinetic energy from fluid to particles
       + (1-killJ2_)*J2
    );

    if (writeTBudget_)
    {
        #include "writeTBudget.H"
    }
    ThetaEqn.relax();
    ThetaEqn.solve();

    // Limit the value of temperature between Tsmall and MaxTheta
    Theta_.max(Tsmall);
    Theta_.min(MaxTheta);

    // Update coefficients after solving Theta
    PsCoeff = granularPressureModel_->granularPressureCoeff(
          alpha_, gs0_, rhoa_, e_);
    pa_ = PsCoeff*Theta_;
    mua_ = viscosityModel_->mua(alpha_+alphaSmall, Theta_+Tsmall, gs0_,
                   muSaltCoef_, K, rhoa_, da_, e_);
    lambda_ = viscosityModel_->lambda(alpha_, Theta_, gs0_, rhoa_, da_, e_);

    pa_.max(0.0);
    pa_.correctBoundaryConditions();

    mua_.max(0.0);
    mua_.correctBoundaryConditions();

    Info << "Granular temp.  Theta: Min ="<< gMin(Theta_)
             << " Max = "<< gMax(Theta_)<< endl;
}

// ************************************************************************* //
