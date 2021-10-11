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
{
    word schemeA("div(alpha,nu)");

    // Compute volume fraction of phase a and b at cell faces
    surfaceScalarField alphaf(fvc::interpolate(alpha));
    surfaceScalarField betaf(scalar(1.0) - alphaf);

    // Compute inverse matrix for momentum balance of phase a and b
    //     at cell centers
    volScalarField rUaA(1.0/UaEqn.A());
    volScalarField rUbA(1.0/UbEqn.A());

    // Interpolate inverse matrix for momentum balance of phase a and b
    //     at cell faces
    surfaceScalarField rUaAf(fvc::interpolate(rUaA));
    surfaceScalarField rUbAf(fvc::interpolate(rUbA));

    // interpolate 1/(alpha+alphaSmall) on the faces
    surfaceScalarField oneOverAlpha
    (
        fvc::interpolate
        (
            scalar(1.0)/(alpha+alphaSmall),
            "interpolate(oneoveralpha)"
        )
    );

    // local Schmidt number
    if (SUSlocal)
    {
        volScalarField ustar
        (
            sqrt(mag(beta*nuEffb*(fvc::grad(Ub)+gradUbT)))
        );
        SUStilde = scalar(1.0) + scalar(2.0)*pow(Ufall/(ustar+Usmall), 2);
        SUStilde.min(scalar(3.0));
        SUS = SUStilde;

        Info<<"ustarmax"<<max(ustar).value()<<endl;
    }

    // PIMPLE algorithm loop
    while (pimple.correct())
    {
        // Compute predicted velocity of phase a and b at cell centers
        volVectorField HabyA(Ua);
        HabyA = rUaA*UaEqn.H();
        volVectorField HbbyA(Ub);
        HbbyA = rUbA*UbEqn.H();

        // Momentum flux at the cell faces for phase a
        surfaceScalarField phiDraga
        (
          // turbulent suspension term
          -(fvc::interpolate
            (
                SUS*(K/rhoa)*turbulenceb->nut()*(beta/(alpha + alphaSmall))
            )
           *fvc::snGrad(alpha))*rUaAf*mesh.magSf()
          // grad(pa)
          - (oneOverAlpha/rhoa*fvc::snGrad(pa))*rUaAf*mesh.magSf()
          // diagonal part of stress --> move to Rca in UEqns.H ?
          + (oneOverAlpha/rhoa*fvc::snGrad(lambdaUa*tr(gradUaT)))
           *rUaAf*mesh.magSf()
          // gravity
          + rUaAf
           *(
               (1.0-rhob/rhoa)*(gravity & mesh.Sf())
            )
          // mean pressure gradient
          + rUaAf*((gradPOSC/(tilt*rhob +(1.0-tilt)*rhoa)) & mesh.Sf())
        );

        // Momentum flux at the cell faces for phase b
        surfaceScalarField phiDragb
        (
          // turbulent suspension term
            fvc::interpolate(SUS*(K/rhob)*turbulenceb->nut())*fvc::snGrad(alpha)
            *rUbAf*mesh.magSf()
          // driving force
          + rUbAf*((gradPOSC/rhob) & mesh.Sf())
          // Convective acceleration for 1D non-uniform simulation
          - rUbAf*fvc::interpolate(neg(alpha-0.3))*((gradConvAcc/rhob)
            & mesh.Sf())
        );

        //    fvOptions.makeRelative(phiHabyA);
        //    fvOptions.makeRelative(phiHbbyA);


        // Fix for gravity on outlet boundary.
        forAll(p_rbgh.boundaryField(), patchi)
        {
             if (isA<zeroGradientFvPatchScalarField>
                 (p_rbgh.boundaryField()[patchi]))
             {
                 phiDraga.boundaryFieldRef()[patchi] = 0.0;
                 phiDragb.boundaryFieldRef()[patchi] = 0.0;
             }
             if (isA<symmetryFvPatchScalarField>
                 (p_rbgh.boundaryField()[patchi]))
             {
                 phiDraga.boundaryFieldRef()[patchi] = 0.0;
                 phiDragb.boundaryFieldRef()[patchi] = 0.0;
             }
             if (p_rbgh.boundaryField().types()[patchi] ==
                 "symmetryPlane")
             {
                 phiDraga.boundaryFieldRef()[patchi] = 0.0;
                 phiDragb.boundaryFieldRef()[patchi] = 0.0;
             }
             if (p_rbgh.boundaryField().types()[patchi] ==
                 "groovyBC")
             {
                 phiDraga.boundaryFieldRef()[patchi] = 0.0;
                 phiDragb.boundaryFieldRef()[patchi] = 0.0;
             }
             if (p_rbgh.boundaryField().types()[patchi] ==
                 "outletInlet")
             {
                 phiDraga.boundaryFieldRef()[patchi] = 0.0;
                 phiDragb.boundaryFieldRef()[patchi] = 0.0;
             }
             if (p_rbgh.boundaryField().types()[patchi] ==
                 "fixedFluxPressure")
             {
                 phiDraga.boundaryFieldRef()[patchi] = 0.0;
                 phiDragb.boundaryFieldRef()[patchi] = 0.0;
             }
        }
        // Compute corrected velocity fluxes at cell faces
        //      (account for time evolution)
        surfaceScalarField phiUa
        (
            (fvc::interpolate(HabyA) & mesh.Sf()) + rUaAf*fvc::ddtCorr(Ua, phia)
        );

        surfaceScalarField phiUb
        (
            (fvc::interpolate(HbbyA) & mesh.Sf()) + rUbAf*fvc::ddtCorr(Ub, phib)
        );
        // Compute predicted velocity fluxes at cell faces
        //      (account for momentum fluxes computed at cell faces)
        surfaceScalarField phiaHbyA(phiUa + phiDraga);
        surfaceScalarField phibHbyA(phiUb + phiDragb);

        // Face-drag coefficients
        surfaceScalarField rAUKda(fvc::interpolate(beta/rhoa*K*rUaA));
        surfaceScalarField rAUKdb(fvc::interpolate(alpha/rhob*K*rUbA));
        // Compute predicted velocity fluxes at cell faces including
        //       including semi-implicit drag contribution
        surfaceScalarField phiHbyA
        (
            alphaf*(phiaHbyA + rAUKda*phib)
          + betaf*(phibHbyA + rAUKdb*phia)
        );

        // adjustPhi(phiHbyA, U, p_rbgh);

        // Assemble Poisson operator for the pressure
        surfaceScalarField Dp
        (
            "(rho*(1|A(U)))",
            alphaf*rUaAf/rhoa + betaf*rUbAf/rhob
        );

        // Update the fixedFluxPressure BCs to ensure flux consistency
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            p_rbgh.boundaryFieldRef(),
            (
                phiHbyA.boundaryField()
              - alphaf.boundaryField()
               *(mesh.Sf().boundaryField() & Ua.boundaryField())
              + betaf.boundaryField()
               *(mesh.Sf().boundaryField() & Ub.boundaryField())
             )
            /(mesh.magSf().boundaryField()*Dp.boundaryField())
        );

        while (pimple.correctNonOrthogonal())
        {
            // Solve for the Poisson equation on the pressure
            fvScalarMatrix pEqn
            (
                fvm::laplacian(Dp, p_rbgh) == fvc::div(phiHbyA)
            );
            // Set pressure reference
            pEqn.setReference(pRefCell, getRefCellValue(p_rbgh, pRefCell));

#if (defined(OFVERSION) && (OFVERSION >= 700))
            pEqn.solve();
#else
            pEqn.solve(mesh.solver(p_rbgh.select(pimple.finalInnerIter())));
#endif

            if (pimple.finalNonOrthogonalIter())
            {
                // Compute pressure gradient at cell faces
                surfaceScalarField SfGradp(pEqn.flux()/Dp);
                // Compute corrected velocity fluxes at cell faces
                //       due to pressure gradient
                surfaceScalarField phias(phiaHbyA - rUaAf*SfGradp/rhoa);
                surfaceScalarField phibs(phibHbyA - rUbAf*SfGradp/rhob);
                // Compute mixture velocity flux at cell faces
                phi = phiHbyA - pEqn.flux();
                // Compute corrected relative velocity fluxes due to drag force
                //       at cell faces
                surfaceScalarField phir
                (
                    ((phias + rAUKda*phibs) - (phibs + rAUKdb*phias))
                   /(1 - rAUKda*rAUKdb)
                );
                // Compute corrected phase velocity fluxes at cell faces
                //       using updated semi-implicit drag
                phia = phi + betaf*(phir);
                phib = phi - alphaf*(phir);

                // Compute contibuity error
                #include "continuityErrs.H"
                p_rbgh.relax();
                // Is it useful?
                SfGradp = pEqn.flux()/Dp;

                // Velocity correction step at cell centers
                volVectorField Ubs
                (
                    HbbyA + rUbA*fvc::reconstruct(phiDragb/rUbAf-SfGradp/rhob)
                );
                volVectorField Uas
                (
                    HabyA + rUaA*fvc::reconstruct(phiDraga/rUaAf-SfGradp/rhoa)
                );
                // Compute corrected relative velocity due to drag force
                //       at cell centers
                volScalarField Da(beta/rhoa*K*rUaA);
                volScalarField Db(alpha/rhob*K*rUbA);
                volVectorField Ur(((1 - Db)*Uas - (1 - Da)*Ubs)/(1 - Da*Db));
                // Compute mixture velocity flux at cell centers
                U = alpha*(Uas + Da*Ub) + beta*(Ubs + Db*Ua);
                // Compute corrected phase velocity at cell centers
                //       using updated semi-implicit drag
                Ua = U + beta*Ur;
                Ua.correctBoundaryConditions();
                fvOptions.correct(Ua);

                Ub = U - alpha*Ur;
                Ub.correctBoundaryConditions();
                fvOptions.correct(Ub);

                // Compute total pressure including hydrostatic contribution
                p = p_rbgh + rhob*gh;
                if (p_rbgh.needReference())
                {
                    p += dimensionedScalar
                    (
                        "p",
                        p.dimensions(),
                        pRefValue - getRefCellValue(p, pRefCell)
                    );
                    p_rbgh = p - rhob*gh;
                    p_rbgh.correctBoundaryConditions();
                }
                // Apply velocity limiters and correct boundary conditions
                #include "UaUbLimiter.H"
                Ua.correctBoundaryConditions();
                // Update mixture velocity
                U = alpha*Ua + beta*Ub;
                U.correctBoundaryConditions();
            }
        }
    }
}
