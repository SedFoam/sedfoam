/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "twophasekOmega.H"

#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(twophasekOmega, 0);
addToRunTimeSelectionTable(RASModel, twophasekOmega, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twophasekOmega::twophasekOmega
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    twophaseRASProperties_
    (
        IOobject
        (
            "twophaseRASProperties",
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    ppProperties_
    (
        IOobject
        (
            "ppProperties",
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    C3om_
    (
        twophaseRASProperties_.lookup("C3om")
    ),
    C4om_
    (
        twophaseRASProperties_.lookup("C4om")
    ),
    KE2_
    (
        twophaseRASProperties_.lookup("KE2")
    ),
    KE4_
    (
        twophaseRASProperties_.lookup("KE4")
    ),    
    omegaBC_
    (
     //twophaseRASProperties_.lookup("omegaBC")
        twophaseRASProperties_.lookupOrDefault("omegaBC",dimensionedScalar("omegaBC",dimensionSet(0,0,-1,0,0,0,0), 0 ))
    ),
    B_
    (
        twophaseRASProperties_.lookup("B")
    ),
    Cmu_
    (
        twophaseRASProperties_.lookup("Cmu")
    ),
    betaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaOmega",
            coeffDict_,
            0.072
        )
    ),
    kSmall_
    (
        twophaseRASProperties_.lookup("kSmall")
    ),
    nutMax_
    (
        twophaseRASProperties_.lookup("nutMax")
    ),
/*    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            coeffDict_,
            0.5
        )
    ),*/
    alphaKOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaKOmega",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            coeffDict_,
            0.52
        )
    ),
    alphaOmegaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmegaOmega",
            coeffDict_,
            0.5
        )
    ),
    alpha_(mesh_.lookupObject<volScalarField> ("alpha")),
    tmfexp_(mesh_.lookupObject<volScalarField> ("tmfexp")),
    ESD3_(mesh_.lookupObject<volScalarField> ("ESD3")),
    ESD4_(mesh_.lookupObject<volScalarField> ("ESD4")),
    ESD5_(mesh_.lookupObject<volScalarField> ("ESD5")),
    ESD_(mesh_.lookupObject<volScalarField> ("ESD")),
    
    
    
    alphaMax_(readScalar(ppProperties_.lookup("alphaMax"))),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    )
{
    bound(k_, kMin_);
    bound(omega_, omegaMin_);

    nut_ = k_/omega_;
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> twophasekOmega::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> twophasekOmega::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}

tmp<fvVectorMatrix> twophasekOmega::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> twophasekOmega::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool twophasekOmega::read()
{
    if (RASModel::read())
    {
        betaOmega_.readIfPresent(coeffDict());
//        alphaK_.readIfPresent(coeffDict());
        alphaKOmega_.readIfPresent(coeffDict());
        alphaOmega_.readIfPresent(coeffDict());
        alphaOmegaOmega_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void twophasekOmega::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    volScalarField G(GName(), nut_*2*magSqr(symm(fvc::grad(U_))));
    //volScalarField G(GName(), fvc::grad(U_) && (-(2.0/3.0)*(k_ + nut_*tr(fvc::grad(U_)))*I + nut_*2*symm(fvc::grad(U_))));

    // syntax for FOAMEXTEND    
//    volScalarField G("RASModel::G", fvc::grad(U_) && (-(2.0/3.0)*(k_ + nut_*tr(fvc::grad(U_)))*I + nut_*2*symm(fvc::grad(U_))));

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::Sp(fvc::div(phi_),omega_)
      - fvm::laplacian(DomegaEff(), omega_, "laplacian(DomegaEff,omega)")
      ==
      - fvm::SuSp (-alphaOmega_*G/max(k_,kSmall_), omega_)
      - fvm::Sp(ESD_,omega_)
      - fvm::Sp(betaOmega_*omega_, omega_)
      + ESD2()*fvm::Sp(C3om_*KE2_, omega_)
      + fvm::Sp((C4om_*KE4_*ESD5_*nut_/k_), omega_)
      // BC in porous bed
      + (-C3om_*KE2_*ESD2() + betaOmega_*omega_ - C4om_*KE4_*ESD5_*nut_/k_)
       *pos(alpha_-0.9*alphaMax_)*omegaBC_
    ); 

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);

    //    Info<<"max1(omega)="<<max(omega_).value()<<endl;

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_, "laplacian(DkEff,k)")
      ==
        //G
      - fvm::SuSp(-G/k_, k_)
      + fvm::Sp(-Cmu_*omega_, k_)
      + fvm::Sp(ESD_, k_)
      + fvm::Sp(KE4_*ESD4_*nut_/k_, k_)
      + ESD2()*fvm::Sp(KE2_, k_)
    );      

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    nut_ = k_/omega_;
    
    nut_.min(nutMax_);
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
