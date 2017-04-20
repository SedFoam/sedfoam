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

surfaceScalarField phiRa =
(
   -fvc::interpolate(nuEffa)*mesh.magSf()*(fvc::snGrad(alpha))
   /fvc::interpolate(alpha+ alphaSmall)
);
tmp<fvVectorMatrix> UaEqn
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
    fvm::Sp(-beta/rhoa*K,Ua)                  // implicit drag term
//+ beta/rhoa*K*Ub                        // explicit drag transfered to p-equation
//- gradPs                                // shear induced pressure gradient moved to p-equation
  - gradPf                                // contact pressure gradient moved to p-equation
  - beta/rhoa*(liftCoeff - Cvm*rhob*DDtUb)  // lift and added mass forces
);
          
UaEqn().relax();

surfaceScalarField phiRb =
(
  - fvc::interpolate(nuEffb)*mesh.magSf()*(fvc::snGrad(beta))
   /fvc::interpolate(beta)
);

tmp<fvVectorMatrix> UbEqn
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
    fvm::Sp(-alpha/rhob*K,Ub)                 // implicit drag term
//+ alpha/rhob*K*Ua                       // Explicit drag transfered to p-equation
  + alpha/rhob*(liftCoeff + Cvm*rhob*DDtUa) // lift and added mass forces
);
UbEqn().relax();
