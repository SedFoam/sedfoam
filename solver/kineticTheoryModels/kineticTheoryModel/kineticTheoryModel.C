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
        kineticTheoryProperties_.lookup("kineticTheory")
    ),
    equilibrium_
    (
        kineticTheoryProperties_.lookup("equilibrium")
    ),
    collisions_
    (
        kineticTheoryProperties_.lookupOrDefault("collisions", false)
    ),
    extended_
    (
        kineticTheoryProperties_.lookupOrDefault("extended", false)
    ),
    softParticles_
    (
        kineticTheoryProperties_.lookupOrDefault("softParticles", false)
    ),
    limitProduction_
    (
        kineticTheoryProperties_.lookupOrDefault("limitProduction", false)
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
    conductivityModelAlpha_
    (
        conductivityModelAlpha::New
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
    frictionalStressModel_
    (
        frictionalStressModel::New
        (
            kineticTheoryProperties_
        )
    ),
    e_
    (
        kineticTheoryProperties_.lookupOrDefault
        (
            "e",
            dimensionedScalar("e",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.9)
        )
    ),
    kn_
    (
        kineticTheoryProperties_.lookupOrDefault
        (
            "kn",
            dimensionedScalar("kn",
                          dimensionSet(1, 0, -2, 0, 0, 0, 0),
                          1e6)
        )
    ),
    muPart_
    (
        kineticTheoryProperties_.lookupOrDefault
        (
            "muPart",
            dimensionedScalar("muPart",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.0)
        )
    ),
    alphaMax_
    (
        kineticTheoryProperties_.lookupOrDefault
        (
            "alphaMax",
            dimensionedScalar("alphaMax",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0.6)
        )
    ),
    DiluteCut_
    (
        kineticTheoryProperties_.lookupOrDefault
        (
            "DiluteCut",
            dimensionedScalar("DiluteCut",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          1e-4)
        )
    ),
    ttzero
    (
        kineticTheoryProperties_.lookupOrDefault
        (
            "ttzero",
            dimensionedScalar("ttzero",
                          dimensionSet(0, 0, 1, 0, 0, 0, 0),
                          0)
        )
    ),
    ttone
    (
        kineticTheoryProperties_.lookupOrDefault
        (
            "ttone",
            dimensionedScalar("ttone",
                          dimensionSet(0, 0, 1, 0, 0, 0, 0),
                          0)
        )
    ),
    MaxTheta
    (
        kineticTheoryProperties_.lookupOrDefault
        (
            "MaxTheta",
            dimensionedScalar("MaxTheta",
                          dimensionSet(0, 2, -2, 0, 0, 0, 0),
                          1e3)
        )
    ),
    phi_
    (
	kineticTheoryProperties_.lookupOrDefault
	(
	    "phi",
	    dimensionedScalar("phi",
			dimensionSet(0, 0, 0, 0, 0, 0, 0),
			32)
	)*M_PI/180.0 //32° angle of repose
    ),
    killJ1_
    (
        kineticTheoryProperties_.lookupOrDefault
        (
            "killJ1_",
            dimensionedScalar
            (
                "killJ1_",
                dimensionSet(0, 0, 0, 0, 0, 0, 0),
                1
            )
        )
    ),
    killFricCorVisc_
    (
        kineticTheoryProperties_.lookupOrDefault("killFricCorVisc", true)
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
    kappa2_
    (
        IOobject
        (
            "kappa2",
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
    if (not kineticTheory_)
    {
        return;
    }

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    const scalar Pi = constant::mathematical::pi;

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

    dimensionedScalar TsmallSqrt = sqrt(Tsmall);

    volScalarField ThetaSqrt = sqrt(Theta_);

    dimensionedScalar Tpsmall_
    (
         "Tpsmall_",
         dimensionSet(1, -1, -3, 0, 0, 0, 0),
         scalar(1e-30)
    );

    // Now start to define terms in the Kinetic theory Equations
    surfaceScalarField phi =
        1.5*rhoa_*phia_*fvc::interpolate((alpha_+alphaSmall));
    surfaceScalarField phicorr = 1.5*alphaSmall*rhoa_*phia_;
//    volVectorField U = fvc::reconstruct(phia_);
    volTensorField dU = fvc::grad(Ua_);//gradUat.T(); //that is fvc::grad(Ua_);
    volSymmTensorField D = symm(dU);   //0.5*(dU + dU.T)
    volSymmTensorField dUU = dev(D);


    // NB, drag = K*alpha*beta,
    // (the alpha and beta has been extracted from the drag function for
    // numerical reasons)
    // this is inconsistent with momentum equation,
    // since the following form missed the drift velocity
    volScalarField Ur_ = mag(Ua_ - Ub_);
    volScalarField betaPrim = (1.0 - alpha_)*draga_.K(Ur_);
    // The following is to calculate parameter tmf_ in u_f*u_s correlation
    volScalarField tmf_ =
        Foam::exp(-B*rhoa_*scalar(6.0)*epsilonb/(max(kb*betaPrim, Tpsmall_)));

    // Calculating the radial distribution function (solid volume fraction is
    //  limited close to the packing limit, but this needs improvements)
    //  The solution is higly unstable close to the packing limit.
    gs0_ = radialModel_->g0
    (
        min(alpha_, alphaMax_ - alphaSmall),
        alphaMax_
    );

    // collision frequency correcion for soft particles 
    volScalarField f_fkt = 1.*(alpha_+alphaSmall)/(alpha_+alphaSmall);

    if (softParticles_)
    {
	//Berzi and Jenkins (2015)
	dimensionedScalar tc = da_/5*pow(rhoa_*Pi*da_/(4*kn_), 1./2.); //Contact duration
	volScalarField fkt = 24/da_*alpha_*gs0_*ThetaSqrt/sqrtPi; //Contact frequency
	f_fkt = 1./(1.+tc*fkt);
    }
    
    // Correction for frcitional particles (Chialvo and Sundaresan (2013))
    dimensionedScalar f_mu = 3./2*muPart_*exp(-3*muPart_);
    dimensionedScalar fric_correction_visc = 1+3./10*pow((1-pow(e_,2)), -2./3)*(1-exp(-8*muPart_));
    dimensionedScalar e_eff = e_ - f_mu;
    dimensionedScalar fric_correction = (1-pow(e_eff,2))/(1-pow(e_, 2));

    // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
    volScalarField PsCoeff = granularPressureModel_->granularPressureCoeff
    (
        alpha_,
        gs0_,
        rhoa_,
        e_
    );
    PsCoeff *= f_fkt; //Correction for soft particles

    volScalarField ThetaClip = Theta_;
    if (limitProduction_)
    {
        // limit the production by clipping Theta to 2/3kb
        volScalarField ThetaClip = min(Theta_, (2.0/3.0)*max(kb));
    }
    // 'thermal' conductivity (Table 3.3, p. 49)
    kappa_ = conductivityModel_->kappa(alpha_, ThetaClip, gs0_, rhoa_, da_, e_);
    kappa_ *= f_fkt; // Correction for soft particles

    // conductivity due to concentration gradient
    kappa2_ = conductivityModelAlpha_->kappa(alpha_+alphaSmall, ThetaClip, gs0_, gs0Prime_, rhoa_, da_, e_);
    kappa2_ *= f_fkt;

    // particle viscosity (Table 3.2, p.47)
    mua_ = viscosityModel_->mua
    (
        alpha_,
        ThetaClip,
        gs0_,
        rhoa_,
        da_,
        e_
    );

    mua_ *= f_fkt; // Correction for soft particles
    if (not killFricCorVisc_){
    mua_ *= fric_correction_visc; // Correction for frictional particles
    }

    // classical kinetic theory
    volScalarField Lc = da_*(alpha_+alphaSmall)/(alpha_+alphaSmall);

    if (extended_)
    {
	volScalarField f2 = rhoa_*da_*
		((4.0/5.0)*sqr(alpha_+alphaSmall)*gs0_*(1.0 + e_)/sqrtPi
      		+ (1.0/15.0)*sqrtPi*gs0_*(1.0 + e_)*(3.0*e_ - 1.0)*sqr(alpha_+alphaSmall)/(3.0 - e_)
      		+ (1.0/6.0)*(alpha_+alphaSmall)*sqrtPi/(3.0 - e_));
	volScalarField f3 = 12./sqrtPi*rhoa_/da_*(1-sqr(e_))*sqr(alpha_+alphaSmall)*gs0_;
	volScalarField Lstar = 1./2*1./2*pow(alpha_*gs0_,1./3)*da_;
        // extended kinetic theory Jenkins (2007)
        Lc = da_*max
        (
            1.,
	    pow(f3/f2,1./3.)*pow(Lstar, 2./3.)
            //0.5*pow(30./(1.+sqr(sqrtPi)/12.)*(1-e_)*sqr(alpha_)*gs0_, 1./3.)
        );
    }


    // bulk viscosity  p. 45 (Lun et al. 1984).
// limit production
    lambda_ = (4.0/3.0)*sqr(alpha_)*rhoa_*da_*gs0_*(1.0+e_)
              *sqrt(ThetaClip)/sqrtPi;

    // stress tensor, Definitions, Table 3.1, p. 43
    volSymmTensorField tau = 2.0*mua_*D + (lambda_ - (2.0/3.0)*mua_)*tr(D)*I;

    if (not equilibrium_)
    {
        // dissipation (Eq. 3.24, p.50)
        volScalarField gammaCoeff =
            (3.0*(1.0 - sqr(e_))*sqr(alpha_)*rhoa_*gs0_
            *((4.0/Lc)*ThetaSqrt/sqrtPi-tr(D)));
	gammaCoeff *= f_fkt; // Correction for soft particles
	gammaCoeff *= fric_correction; // Correction for frictional particles
        // Eq. 3.25, p. 50 Js = J1 - J2
        volScalarField J1 = 1.8*3.0*alpha_*betaPrim;
        /*
        volScalarField J2 =
            0.25*alpha_*sqr(betaPrim)*da_*sqr(Ur_)
           /(rhoa_*sqrtPi*ThetaSqrt);
        */


        // construct the granular temperature equation (Eq. 3.20, p. 44)
        // NB. note that there are two typos in Eq. 3.20
        // no grad infront of Ps
        // wrong sign infront of laplacian
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
           // granular temperature conduction due to concentration gradient.
          + fvc::laplacian(kappa2_, alpha_, "laplacian(kappa2,alpha)")
           // energy disipation due to inelastic collision.
          + fvm::Sp(-gammaCoeff, Theta_)
           // dissipation due to interphase slip
          + fvm::Sp(-J1, Theta_)
           // turbulence generation due to fluid turbulence
          + (1-killJ1_)*scalar(2.0/3.0)*J1*tmf_*kb
        );

        ThetaEqn.relax();
        ThetaEqn.solve();
    }
    else
    {
        // equilibrium => dissipation == production
        // Eq. 4.14, p.82
        volScalarField K1 = 2.0*(1.0 + e_)*rhoa_*gs0_;
        volScalarField K3 = 0.5*da_*rhoa_*
            (
                (sqrtPi/(3.0*(3.0-e_)))
               *(1.0 + 0.4*(1.0 + e_)*(3.0*e_ - 1.0)*alpha_*gs0_)
               +1.6*alpha_*gs0_*(1.0 + e_)/sqrtPi
            );

        volScalarField K2 =
            4.0*da_*rhoa_*(1.0 + e_)*alpha_*gs0_/(3.0*sqrtPi) - 2.0*K3/3.0;

        volScalarField K4 = 12.0*(1.0 - sqr(e_))*rhoa_*gs0_/(da_*sqrtPi);

        volScalarField trD = tr(D);
        volScalarField tr2D = sqr(trD);
        volScalarField trD2 = tr(D & D);

        volScalarField t1 = K1*alpha_;
        volScalarField l1 = -t1*trD;
        volScalarField l2 = sqr(t1)*tr2D;
        volScalarField l3 = 4.0*K4*alpha_*(2.0*K3*trD2 + K2*tr2D);

        Theta_ = sqr((l1 + sqrt(l2 + l3))/(2.0*(alpha_ + 1.0e-4)*K4));
    }

    // Clipping of Theta for stability by Z. Cheng
    //    (need to check if necessary...)
    forAll(alpha_, cellk)
    {
// for initial stability
       if (tt.value()<=ttzero.value() && alpha_[cellk]>= 0)
       {
            Theta_[cellk] = 0.0*Theta_[cellk];
       }
       if (tt.value()<=ttone.value() && alpha_[cellk] >= 0.5)
       {
            Theta_[cellk] = 0.0*Theta_[cellk];
       }
       if (tt.value()>ttone.value() && alpha_[cellk]<=alphaSmall.value())
       {
            Theta_[cellk] =
                min(Theta_[cellk], scalar(2.0)/scalar(3.0)*kb[cellk]);
       }
    }

    Theta_.max(Tsmall);
    Theta_.min(MaxTheta);

// need to update after solving Theta Equation.
    PsCoeff = granularPressureModel_->granularPressureCoeff
    (
        alpha_,
        gs0_,
        rhoa_,
        e_
    );
    // update bulk viscosity and shear viscosity
    // bulk viscosity  p. 45 (Lun et al. 1984).
    lambda_ = (4.0/3.0)*sqr(alpha_)*rhoa_*da_*gs0_*(1.0+e_)*sqrt(Theta_)/sqrtPi;
    // particle viscosity (Table 3.2, p.47)
    mua_ = viscosityModel_->mua(alpha_, Theta_, gs0_, rhoa_, da_, e_);


    // update particle pressure, just from the kinetic part
    pa_ = PsCoeff*Theta_;

    mua_.max(0.0);
    mua_.correctBoundaryConditions();

    Info << "Granular temp.  Theta: Min ="<< gMin(Theta_)
             << " Max = "<< gMax(Theta_)<< endl;
}
// test function for callFrictionStress.H
void Foam::kineticTheoryModel::updateRheo
(
    const volScalarField& kb,
    const volScalarField& epsilonb,
    const dimensionedScalar& B
)
{
    const scalar sqrtPi = sqrt(constant::mathematical::pi);

    dimensionedScalar alphaSmall
    (
        "small",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        scalar(1.0e-4)
    );

    dimensionedScalar Tpsmall_
    (
        "Tpsmall_",
        dimensionSet(1, -1, -3, 0, 0, 0, 0),
        scalar(1e-30)
    );

    dimensionedScalar Kasmall_
    (
        "Kasmall_",
        dimensionSet(1, -3, -1, 0, 0, 0, 0),
        scalar(1e-12)
    );

     // NB, drag = K*alpha*beta,
    // (the alpha and beta has been extracted from the drag function for
    // numerical reasons)
    // this is inconsistent with momentum equation,
    // since the following form missed the drift velocity
    volScalarField Ur_ = mag(Ua_ - Ub_);
    volScalarField Ka = draga_.K(Ur_);
    volScalarField betaPrim = (1.0 - alpha_)*Ka;
    // The following is to calculate parameter tmf_ in u_f*u_s correlation
    volScalarField tmf_ =
        Foam::exp(-B*rhoa_*scalar(6.0)*epsilonb/(max(kb*betaPrim, Tpsmall_)));

    gs0_ = radialModel_->g0
    (
        min(alpha_, alphaMax_ - alphaSmall),
        alphaMax_
    );

    volScalarField F_ = (1.0/6.0)*alpha_*sqrtPi;

    //volScalarField F_ = 0.12*alpha_*Foam::exp(-600*pow(alpha_,5));
    //volScalarField F_ = 0.12*alpha_*Foam::exp(-600*pow(alpha_,7));


    if (collisions_)
    {
        volScalarField dum =
            1.0/(1.0+((1.0-e_*e_)/(Ka+Kasmall_))
            *(4*rhoa_*alpha_*gs0_/(da_*sqrtPi))*pow(Theta_, 0.5));
        Theta_ = dum*(2.0/3.0)*tmf_*kb;
        Info << "run avec collisions " << endl;
    }

    else
    {
        Theta_ = (2.0/3.0)*tmf_*kb;
    }

    mua_= rhoa_*F_*da_*pow(Theta_, 0.5);
    
    pa_ = rhoa_*alpha_*Theta_;
}
//}

// ************************************************************************* //
