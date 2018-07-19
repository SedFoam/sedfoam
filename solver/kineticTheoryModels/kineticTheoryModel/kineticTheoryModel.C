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
    kineticTheory_(kineticTheoryProperties_.lookup("kineticTheory")),
    equilibrium_(kineticTheoryProperties_.lookup("equilibrium")),
    collisions_(kineticTheoryProperties_.lookupOrDefault("collisions_", false)),
    extended_(kineticTheoryProperties_.lookupOrDefault("extended_", false)),
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
    e_(kineticTheoryProperties_.lookup("e")),
    alphaMax_(kineticTheoryProperties_.lookup("alphaMax")),
    DiluteCut_(kineticTheoryProperties_.lookup("DiluteCut")),
    ttzero(kineticTheoryProperties_.lookup("ttzero")),
    ttone(kineticTheoryProperties_.lookup("ttone")),
    MaxTheta(kineticTheoryProperties_.lookup("MaxTheta")),
    phi_(dimensionedScalar(kineticTheoryProperties_.lookup("phi"))*M_PI/180.0),
    //    relaxPaKin_(kineticTheoryProperties_.lookupOrDefault(("relaxPaKin"),1.0),
    relaxPaKin_(kineticTheoryProperties_.lookupOrDefault("relaxPaKin_",dimensionedScalar("relaxPaKin_",dimensionSet(0,0,0,0,0,0,0), 1 ))),
    Theta_
    (
        IOobject
        (
            "Theta",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        Ua_.mesh()
    ),
    mua_
    (
        IOobject
        (
            "mua",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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
            IOobject::AUTO_WRITE
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

void Foam::kineticTheoryModel::solve(const surfaceScalarField& galpha, const volTensorField& gradUat, const volScalarField& kb,const volScalarField& epsilonb,const volScalarField& nutf,const dimensionedScalar& B, const dimensionedScalar& tt)
{
    if (!kineticTheory_)
    {
        return;
    }

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    dimensionedScalar alphaSmall
    (
        "small",
        dimensionSet(0 , 0 ,0 ,0 , 0, 0, 0),
        scalar(1.0e-4)
    );
    dimensionedScalar Tsmall
    (
        "small",
        dimensionSet(0 , 2 ,-2 ,0 , 0, 0, 0),
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
    surfaceScalarField phi = 1.5*rhoa_*phia_*fvc::interpolate((alpha_+alphaSmall));
    surfaceScalarField phicorr = 1.5*alphaSmall*rhoa_*phia_;
//    volVectorField U = fvc::reconstruct(phia_);
    volTensorField dU = fvc::grad(Ua_);//gradUat.T();         //that is fvc::grad(Ua_);
    volSymmTensorField D = symm(dU);         //0.5*(dU + dU.T)
    volSymmTensorField dUU = dev(D);


    // NB, drag = K*alpha*beta,
    // (the alpha and beta has been extracted from the drag function for
    // numerical reasons)
    // this is inconsistent with momentum equation, since the following form missed
    // the drift velocity
    volScalarField Ur_ = mag(Ua_ - Ub_);
    volScalarField betaPrim = (1.0 - alpha_)*draga_.K(Ur_);
    // The following is to calculate parameter tmf_ in u_f*u_s correlation
    volScalarField tmf_ = Foam::exp(-B*rhoa_*scalar(6.0)*epsilonb/(max(kb*betaPrim,Tpsmall_)));

    // Calculating the radial distribution function (solid volume fraction is
    //  limited close to the packing limit, but this needs improvements)
    //  The solution is higly unstable close to the packing limit.
    gs0_ = radialModel_->g0
    (
        min(alpha_, alphaMax_ - alphaSmall),
        alphaMax_
    );

    // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
    volScalarField PsCoeff = granularPressureModel_->granularPressureCoeff
    (
        alpha_,
        gs0_,
        rhoa_,
        e_
    );
 
    // 'thermal' conductivity (Table 3.3, p. 49)
    kappa_ = conductivityModel_->kappa(alpha_, Theta_, gs0_, rhoa_, da_, e_);

    // particle viscosity (Table 3.2, p.47)
// limit the production
    mua_ = viscosityModel_->mua(alpha_,min(Theta_,(2.0/3.0)*max(kb)), gs0_, rhoa_, da_, e_);    

    // classical kinetic theory
    volScalarField Lc = da_*(alpha_+alphaSmall)/(alpha_+alphaSmall);

    if (extended_)
    {
        // extended kinetic theory Jenkins (2007)
        Lc = da_*max(1.,0.5*pow(30./(1.+sqr(sqrtPi)/12.)*(1-e_)*sqr(alpha_)*gs0_,1./3.));
    }

    // dissipation (Eq. 3.24, p.50)
    volScalarField gammaCoeff =
                3.0*(1.0 - sqr(e_))*sqr(alpha_)*rhoa_*gs0_*((4.0/Lc)*ThetaSqrt/sqrtPi-tr(D));

    // Eq. 3.25, p. 50 Js = J1 - J2
    volScalarField J1 = 3.0*alpha_*betaPrim;
    volScalarField J2 =
        0.25*alpha_*sqr(betaPrim)*da_*sqr(Ur_)
       /(rhoa_*sqrtPi*ThetaSqrt);

    // bulk viscosity  p. 45 (Lun et al. 1984).
// limit production
    lambda_ = (4.0/3.0)*sqr(alpha_)*rhoa_*da_*gs0_*(1.0+e_)*sqrt(min(Theta_,(2.0/3.0)*max(kb)))/sqrtPi;    
    // stress tensor, Definitions, Table 3.1, p. 43
    volSymmTensorField tau = 2.0*mua_*D + (lambda_ - (2.0/3.0)*mua_)*tr(D)*I;

    if (!equilibrium_)
    {
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
          fvm::SuSp(-PsCoeff*fvc::div(phia_),Theta_)
           // production due to shear.
          + (tau && dU)
           // granular temperature conduction.
          + fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
           // energy disipation due to inelastic collision.
          + fvm::Sp(-gammaCoeff, Theta_)
           // dissipation due to interphase slip
          + fvm::Sp(-J1, Theta_)
           // turbulence generation due to fluid turbulence
          + scalar(2.0/3.0)*J1*tmf_*kb
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
          Theta_[cellk] = min(Theta_[cellk],scalar(2.0)/scalar(3.0)*kb[cellk]);
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

    /* pf_ is now computed outside the  1.19877
kineticTheory module, see callFrictionStress.H 
    pf_ = frictionalStressModel_->frictionalPressure
    (
        alpha_,
        alphaMinFriction_,
        alphaMax_,
        Fr_,
        eta_,
        p_
    );
    */
//  yes, after solving Theta, because frictional stress is not a part of kinetic theoty
//    PsCoeff = PsCoeff + pf/(max(Theta_,Tsmall));

//    PsCoeff.min(1e+6);
//    PsCoeff.max(-1e+6);

    // update particle pressure, just from the kinetic part
    pa_ = PsCoeff*Theta_;


    /* muf_ is now computed outside the module, see callFrictionStress.H
    // frictional shear stress, Eq. 3.30, p. 52
    muf_ = frictionalStressModel_->muf
    (
        alpha_,
        Theta_,
        alphaMinFriction_, 1.19877

        alphaMax_,
        pf_,
        D,
        phi_
    );
    */

   // add frictional stress
    //mua_ = mua_+muf;
    mua_.max(0.0);
    mua_.correctBoundaryConditions();
    //    muf_.max(0.0); // computed in callFrictionStress.H

    Info<< "kinTheory: max(Theta) = " << max(Theta_).value() <<" min(Theta) = "<<min(Theta_).value()<< endl;

    /* ppMagf is computed in callFrictionStress.H 
    volScalarField ppMagfFriction = frictionalStressModel_->frictionalPressurePrime
    (
        alpha_,
        alphaMinFriction_,
        alphaMax_,
        Fr_,
        eta_,
        p_
    );
    forAll(alpha_, cellI)
    {
        if (alpha_[cellI] >= alphaMinFriction_.value())
        {
            ppMagf_[cellI] = ppMagfFriction[cellI];
        }
    }
    ppMagf_.max(0.0);
    */

}

// test function for callFrictionStress.H
void Foam::kineticTheoryModel::updateRheo (const volScalarField& kb,const volScalarField& epsilonb, const dimensionedScalar& B)
{  
    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    

    dimensionedScalar alphaSmall
    	(
       	   "small",
           dimensionSet(0 , 0 ,0 ,0 , 0, 0, 0),
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
    // this is inconsistent with momentum equation, since the following form missed
    // the drift velocity
    volScalarField Ur_ = mag(Ua_ - Ub_);
    volScalarField Ka = draga_.K(Ur_);
    volScalarField betaPrim = (1.0 - alpha_)*Ka;
    // The following is to calculate parameter tmf_ in u_f*u_s correlation
    volScalarField tmf_ = Foam::exp(-B*rhoa_*scalar(6.0)*epsilonb/(max(kb*betaPrim,Tpsmall_)));
    
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
     volScalarField dum = 1.0/(1.0+((1.0-e_*e_)/(Ka+Kasmall_))*(4*rhoa_*alpha_*gs0_/(da_*sqrtPi))*pow(Theta_,0.5)); 
     Theta_ = dum*(2.0/3.0)*tmf_*kb;
     Info<< "run avec collisions " << endl;
   }
   
   else
   {
   Theta_ = (2.0/3.0)*tmf_*kb;
   }

   mua_= rhoa_*F_*da_*pow(Theta_,0.5);
   //volScalarField paOld = pa_;
   pa_ = rhoa_*alpha_*Theta_;
   //pa_ = paOld + relaxPaKin_*(pa_ - paOld);
   //Info<< "Theta = " << Theta_ << endl;
   //Info<< "mua = " << mua_ << endl;
}




   
//} 

// ************************************************************************* //
