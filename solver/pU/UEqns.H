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
scalar pi(constant::mathematical::pi);

// following is general form of mean flow + sine wave and 2nd order stokes waves
// Add oscilitary pressure gradient
gradPOSC =
(
    gradPMEAN
  + gradPAMP1*Foam::cos(initTheta1*(pi)/180.0 +  2.0*(pi)*runTime.value()/oscpT)
  + gradPAMP2*Foam::cos(initTheta2*(pi)/180.0 +  4.0*(pi)*runTime.value()/oscpT)
  + gradPAMP3*Foam::cos(initTheta3*(pi)/180.0 +  6.0*(pi)*runTime.value()/oscpT)
  + gradPAMP4*Foam::cos(initTheta4*(pi)/180.0 +  8.0*(pi)*runTime.value()/oscpT)
  + gradPAMP5*Foam::cos(initTheta5*(pi)/180.0 + 10.0*(pi)*runTime.value()/oscpT)
);

// Compute mixture density
rho = alpha*rhoa+beta*rhob;

// Compute concentration gradient stress contribution on cell faces for phase a
surfaceScalarField phiRa
(
   -fvc::interpolate(nuEffa)*mesh.magSf()*(fvc::snGrad(alpha))
   /fvc::interpolate(alpha+ alphaSmall)
);

// Assemble momentum balance equation for phase a
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
  + fvm::div(phiRa, Ua, "div(phiRa,Ua)")
  - fvm::Sp(fvc::div(phiRa), Ua)
    ==
  - beta/rhoa*(liftCoeff - Cvm*rhob*DDtUb)  // lift and added mass forces
  + fvOptions(alpha, Ua)
);
UaEqn.relax();
UaEqn += fvm::Sp(beta/rhoa*K, Ua);
fvOptions.constrain(UaEqn);
Ua.correctBoundaryConditions();
fvOptions.correct(Ua);

// Compute concentration gradient stress contribution on cell faces for phase b
surfaceScalarField phiRb
(
  - fvc::interpolate(nuEffb)*mesh.magSf()*(fvc::snGrad(beta))
   /fvc::interpolate(beta)
);

// Assemble momentum balance equation for phase b
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
  + fvOptions(beta, Ub)
);
UbEqn.relax();
UbEqn += fvm::Sp(alpha/rhob*K, Ub);
fvOptions.constrain(UbEqn);
Ub.correctBoundaryConditions();
fvOptions.correct(Ub);
