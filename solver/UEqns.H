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
scalar pi = constant::mathematical::pi;

// following is general form of mean flow + sine wave and 2nd order stokes waves
// Add oscilitary pressure gradient
gradPOSC =
(
    gradPMEAN
  + gradPAMP1*Foam::cos(initTheta*(pi)/180.0 + 2.0*(pi)*runTime.value()/oscpT)
  + gradPAMP2*Foam::sin(initTheta*(pi)/90.0 + 4.0*(pi)*runTime.value()/oscpT)
);

rho = alpha*rhoa+beta*rhob;


/* JC Contribution from Hugo commented
//added tensor (comment if not added)
dimensionedTensor iso
(
  "iso",
  dimensionSet(0, 0, 0, 0, 0, 0, 0),
  tensor(1e0, 0e0 , 0e0, 0e0, 1e0, 0e0, 0e0, 0e0 ,1e0)
  //      XX, XY,   XZ,   YX, YY,  YZ,  ZX,  ZY,  ZZ
);

dimensionedTensor aniso
(
  "aniso",
  dimensionSet(0, 0, 0, 0, 0, 0, 0),
  tensor(0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0 ,0.0)
  //      XX, XY,  XZ,  YX,  YY,  YZ,  ZX,  ZY,  ZZ
);

volVectorField susa
(
    "susa",
    //turbulence suspension
    SUS*(K/rhoa)*beta/(alpha+alphaSmall)*turbulenceb->nut()
    *((SUS_I*iso-SUS_A*aniso) & fvc::grad(alpha))
);
volVectorField susb
(
    "susb",
    //turbulence suspension
    SUS*(K/rhob)*turbulenceb->nut()*((SUS_I*iso-SUS_A*aniso) & fvc::grad(alpha))
);
//end added tensor
// end contribution from Hugo
    */


surfaceScalarField phiRa =
(
   -fvc::interpolate(nuEffa)*mesh.magSf()*(fvc::snGrad(alpha))
   /fvc::interpolate(alpha+ alphaSmall)
);




fvVectorMatrix UaEqn
(
    (scalar(1.0) + Cvm*rhob*beta/rhoa)*
    (
        fvm::ddt(Ua)
      + fvm::div(phia, Ua, "div(phia,Ua)")
      - fvm::Sp(fvc::div(phia), Ua)
    )
// effective stress term from dispersed phase
  - scalar(1.0)/(alpha+alphaSmall)*fvm::laplacian(nuFra, Ua)
//- fvm::laplacian(nuFra, Ua)
  - fvm::laplacian(nuEffa, Ua)
  + fvm::div(phiRa, Ua, "div(phiRa,Ua)")
  - fvm::Sp(fvc::div(phiRa), Ua)
    ==
//    fvm::Sp(-beta/rhoa*K,Ua)        // implicit drag term
//+ beta/rhoa*K*Ub                    // explicit drag transfered to p-equation
    // Hugo  -susa                    // term to add with tensor
//- gradPs              // shear induced pressure gradient moved to p-equation
//  - gradPf                  // contact pressure gradient moved to p-equation
  - beta/rhoa*(liftCoeff - Cvm*rhob*DDtUb)  // lift and added mass forces
  //    + fvOptions(Ua)
);
UaEqn.relax();
UaEqn += fvm::Sp(beta/rhoa*K, Ua);
//fvOptions.constrain(UaEqn());

surfaceScalarField phiRb =
(
  - fvc::interpolate(nuEffb)*mesh.magSf()*(fvc::snGrad(beta))
   /fvc::interpolate(beta)
);

fvVectorMatrix UbEqn
(
    (scalar(1.0) + Cvm*rhob*alpha/rhob)*
    (
        fvm::ddt(Ub)
      + fvm::div(phib, Ub, "div(phib,Ub)")
      - fvm::Sp(fvc::div(phib), Ub)
    )
  - fvm::laplacian(nuEffb, Ub)
  + fvm::div(phiRb, Ub, "div(phiRb,Ub)")
  - fvm::Sp(fvc::div(phiRb), Ub)
    ==
//    fvm::Sp(-alpha/rhob*K,Ub)                 // implicit drag term
    // Hugo   susb                                 //term to add with tensor
//+ alpha/rhob*K*Ua                 // Explicit drag transfered to p-equation
    alpha/rhob*(liftCoeff + Cvm*rhob*DDtUa) // lift and added mass forces
//  + fvOptions(Ub)
);
UbEqn.relax();
UbEqn += fvm::Sp(alpha/rhob*K, Ub);

//fvOptions.constrain(UbEqn());

/*    if (pimple.momentumPredictor())
    {
    word schemeA("div(alpha,nu)");

// limiter for the explicit stress terms to avoid numerical instability and/or
// too strong limitation on the time step
    volScalarField nuFraExp = nuFra;
    nuFraExp.min(nub*1e3);
    nuFraExp.correctBoundaryConditions();

    surfaceScalarField oneOverAlpha =
    (
        fvc::interpolate
        (
            scalar(1.0)/(alpha+alphaSmall),
            "interpolate(oneoveralpha)"
        )
    );


        solve
        (
            UbEqn()
         ==
            alpha/rhob*K*Ua
          - fvc::grad(p_rbgh)/rhob
          + (SUS*(K/rhob)*turbulenceb->nut())*fvc::grad(alpha)
          - (scalar(2.0)/(beta*scalar(3.0)))
            *fvc::grad(((nuEffb*tr(gradUbT)))*beta)
          + (scalar(1.0)/beta)
            *(fvc::div(beta*nuEffb*gradUbT,schemeA))
          + (gradPOSC/rhob)
        );
        solve
        (
            UaEqn()
         ==
            beta/rhoa*K*Ub
          - fvc::grad(p_rbgh)/rhoa
          fvc::reconstruct
          (
          )
        );
    }
*/
