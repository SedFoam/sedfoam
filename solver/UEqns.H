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
// Compute mixture density
rho = alpha*rhoa+beta*rhob;

// Solid phase momentum balance equation
volTensorField Rca
(
    "Rca",
    (((2.0/3.0)*nuEffa-lambdaUa/rhoa)*I*tr(gradUaT) 
    + ((2.0/3.0)*I)*sqr(Ct)*turbulenceb->k() 
    - nuEffa*gradUaT
    )
);

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
  - fvm::laplacian(nuEffa, Ua)
  - fvm::Sp(fvc::div(phiRa), Ua)
  + fvm::div(phiRa, Ua, "div(phiRa,Ua)")
  + fvc::div(Rca)
  + (fvc::grad(alpha)/(alpha+alphaSmall) & Rca)
    ==
  - beta/rhoa*(liftCoeff - Cvm*rhob*DDtUb)  // lift and added mass forces
  //    + fvOptions(alpha,Ua)
);
UaEqn.relax();
UaEqn += fvm::Sp(beta/rhoa*K, Ua);
//fvOptions.constrain(UaEqn());
//Ua.correctBoundaryConditions();
//fvOptions.correct(Ua);

// Fluid phase momentum balance equation
volTensorField Rcb
(
    "Rcb",
    ((2.0/3.0)*I)*nuEffb*tr(gradUbT) - nuEffb*gradUbT
);
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
    alpha/rhob*(liftCoeff + Cvm*rhob*DDtUa) // lift and added mass forces
//  + fvOptions(Ub)
);
UbEqn.relax();
UbEqn += fvm::Sp(alpha/rhob*K, Ub);

//fvOptions.constrain(UbEqn());
//Ub.correctBoundaryConditions();
//fvOptions.correct(Ub);
