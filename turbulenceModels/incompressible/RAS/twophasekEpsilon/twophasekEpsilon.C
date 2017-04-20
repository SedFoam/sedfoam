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

#include "twophasekEpsilon.H"
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

defineTypeNameAndDebug(twophasekEpsilon, 0);
addToRunTimeSelectionTable(RASModel, twophasekEpsilon, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twophasekEpsilon::twophasekEpsilon
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
        C3ep_
    (
        twophaseRASProperties_.lookup("C3ep")
    ),
    C4ep_
    (
        twophaseRASProperties_.lookup("C4ep")
    ),
    KE2_
    (
        twophaseRASProperties_.lookup("KE2")
    ),
    KE4_
    (
        twophaseRASProperties_.lookup("KE4")
    ),    
    B_
    (
        twophaseRASProperties_.lookup("B")
    ),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
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
    alphak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphak",
            coeffDict_,
            1.0
        )
    ),
        alphaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaEps",
            coeffDict_,
	    1.3
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
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateEpsilon("epsilon", mesh_)
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
    bound(epsilon_, epsilonMin_);

    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> twophasekEpsilon::R() const
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


tmp<volSymmTensorField> twophasekEpsilon::devReff() const
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


tmp<fvVectorMatrix> twophasekEpsilon::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> twophasekEpsilon::divDevRhoReff
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


bool twophasekEpsilon::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
	// sigmaEps_.readIfPresent(coeffDict());
	alphak_.readIfPresent(coeffDict());
        alphaEps_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void twophasekEpsilon::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    //volScalarField G(GName(), nut_*2*magSqr(symm(fvc::grad(U_))));
    volScalarField G(GName(), fvc::grad(U_) && (-(2.0/3.0)*(k_ + nut_*tr(fvc::grad(U_)))*I + nut_*2*symm(fvc::grad(U_))));

    forAll(G.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(G.boundaryField()[patchi]))
        {
            G.boundaryField()[patchi] = 0.0;
        }
    }


    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();


    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::Sp(fvc::div(phi_),epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
      - fvm::SuSp(-C1_*G/k_, epsilon_)
      + fvm::Sp(-C2_*epsilon_/k_, epsilon_)
      + fvm::Sp(C3ep_*ESD_,epsilon_)
      + ESD2()*fvm::Sp(C3ep_*KE2_,epsilon_)
      + fvm::Sp(C4ep_*KE4_*ESD5_*nut_/k_,epsilon_)
      //+ fvm::Sp((C4ep_*KE4_*ESD5_*nut_/max(k_,kSmall_)),epsilon_)
    );

    epsEqn().relax();

    epsEqn().boundaryManipulate(epsilon_.boundaryField());

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
      - fvm::SuSp(-G/k_, k_)
      - fvm::Sp(epsilon_/k_, k_)
      + fvm::Sp(ESD_, k_)
      + fvm::Sp(KE4_*ESD4_*nut_/k_, k_)
      + ESD2()*fvm::Sp(KE2_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.min(nutMax_);
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
