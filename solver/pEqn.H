{
    surfaceScalarField alphaf = fvc::interpolate(alpha);
    surfaceScalarField betaf = scalar(1.0) - alphaf;

    volScalarField rUaA = 1.0/UaEqn().A();
    volScalarField rUbA = 1.0/UbEqn().A();

    rUaAf = fvc::interpolate(rUaA);
    surfaceScalarField rUbAf = fvc::interpolate(rUbA);

    volVectorField HabyA = Ua;
    HabyA = rUaA*UaEqn().H();
    volVectorField HbbyA = Ub;
    HbbyA = rUbA*UbEqn().H();
    word schemeA("div(alpha,nu)");
 
    // limiter for the explicit stress terms to avoid numerical instability and/or
    // too strong limitation on the time step
    volScalarField nuFraExp = nuFra;
    nuFraExp.min(nub*1e3);
    nuFraExp.correctBoundaryConditions();

    // interpolate 1/(alpha+alphaSmall) on the faces
    surfaceScalarField oneOverAlpha =
    (
        fvc::interpolate
        (
            scalar(1.0)/(alpha+alphaSmall),
            "interpolate(oneoveralpha)"
        )
    );

    // Momentum fluxes computed on the faces for the particle phase phiDraga and
    // the fluid phase phiDragb
    surfaceScalarField phiDraga =
    (
    // explicit drag term
        fvc::interpolate(beta/rhoa*K*rUaA)*phib 
    // turbulent suspension term
      - (fvc::interpolate(SUS*(K/rhoa)*turbulenceb->nut()*(beta/(alpha + alphaSmall)))
       *fvc::snGrad(alpha))*rUaAf*mesh.magSf() 
    // grad(pa)
      - (oneOverAlpha/rhoa*fvc::snGrad(pa))*rUaAf*mesh.magSf()
    // grad(Pf)
    // - (oneOverAlpha/rhoa*fvc::snGrad(pff))*rUaAf*mesh.magSf() 
    // diagonal part of stress
      + (oneOverAlpha/rhoa*fvc::snGrad(lambdaUa*tr(gradUaT)))*rUaAf*mesh.magSf()
    // diagonal part of stress 2
      - (
            scalar(2.0/3.0)*oneOverAlpha
           *fvc::snGrad((alpha+alphaSmall)*(sqr(Ct)*turbulenceb->k()
          + nuEffa*tr(gradUaT)))
        )
       *rUaAf*mesh.magSf() 
      - (scalar(2.0/3.0)*oneOverAlpha*fvc::snGrad(nuFraExp*tr(gradUaT)))
       *rUaAf*mesh.magSf() 
    // stress part 3
      + rUaAf*(oneOverAlpha
       *fvc::interpolate(fvc::div((alpha+alphaSmall)*nuEffa*gradUaT, schemeA))
      & mesh.Sf())
      + rUaAf*(oneOverAlpha
       *fvc::interpolate(fvc::div(nuFraExp*gradUaT, schemeA)) & mesh.Sf())
    // gravity
      + rUaAf*(g & mesh.Sf())
    // mean pressure gradient
      + rUaAf*((gradPOSC/(tilt*rhob +(1.0-tilt)*rhoa)) & mesh.Sf())
    );
       
    surfaceScalarField phiDragb =
    (
// explicit drag term
        fvc::interpolate(alpha/rhob*K*rUbA)*phia 
// turbulent suspension term
      + (fvc::interpolate(SUS*(K/rhob)*turbulenceb->nut())*fvc::snGrad(alpha))
       *rUbAf*mesh.magSf() 
// diagnoal part of stress
      - (fvc::interpolate(scalar(2.0)/(beta*scalar(3.0)))
       *fvc::snGrad(((nuEffb*tr(gradUbT)))*beta))*rUbAf*mesh.magSf() 
// stress part 3 
      + rUbAf*(fvc::interpolate(scalar(1.0)/beta)
       *fvc::interpolate(fvc::div(beta*nuEffb*gradUbT,schemeA)) & mesh.Sf())
// Gravity force
       +rUbAf*(g & mesh.Sf())
// driving force
      + rUbAf*((gradPOSC/rhob) & mesh.Sf())
    );


    // Fix for gravity on outlet boundary.
    forAll(p.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
        if (isA<symmetryFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
        if (p.boundaryField().types()[patchi] == "symmetryPlane")
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
        if (p.boundaryField().types()[patchi] == "groovyBC")
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
	if (p.boundaryField().types()[patchi] == "outletInlet")
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
    }

    surfaceScalarField phiUa =
    (
        (fvc::interpolate(HabyA) & mesh.Sf()) + rUaAf*fvc::ddtCorr(Ua, phia)
    );

    surfaceScalarField phiUb =
    (
        (fvc::interpolate(HbbyA) & mesh.Sf()) + rUbAf*fvc::ddtCorr(Ub, phib)
    );

    surfaceScalarField phiaHbyA = phiUa + phiDraga;
    surfaceScalarField phibHbyA = phiUb + phiDragb;

    surfaceScalarField phiHbyA = alphaf*phiaHbyA + betaf*phibHbyA;

    adjustPhi(phiHbyA, U, p);

    // Poisson equation for the pressure 
    surfaceScalarField Dp
    (
        "(rho*(1|A(U)))",
        alphaf*rUaAf/rhoa + betaf*rUbAf/rhob
    );

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p.boundaryField(),
        (
            phiHbyA.boundaryField()
          - (alphaf.boundaryField()*(mesh.Sf().boundaryField() & Ua.boundaryField())
          + betaf.boundaryField()*(mesh.Sf().boundaryField() & Ub.boundaryField())
            )
         )
        /(mesh.magSf().boundaryField()*Dp.boundaryField())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField SfGradp = pEqn.flux()/Dp;

            phia = phiaHbyA - rUaAf*SfGradp/rhoa;
            phib = phibHbyA - rUbAf*SfGradp/rhob;
            phi = phiHbyA - pEqn.flux();
            #include "continuityErrs.H"
            // Velocity correction step    
            Ub = HbbyA + rUbA*fvc::reconstruct(phiDragb/rUbAf-SfGradp/rhob);
            Ub.correctBoundaryConditions();
            Ua = HabyA + rUaA*fvc::reconstruct(phiDraga/rUaAf-SfGradp/rhoa);
            #include "UaUbLimiter.H"
            Ua.correctBoundaryConditions();
            // Update mixture velocity    
            U = alpha*Ua + beta*Ub;
            U.correctBoundaryConditions();
        }
    }
}
