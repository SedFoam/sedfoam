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
    word scheme("div(phi,alpha)");
    word schemer("div(phir,alpha)");

    surfaceScalarField phic("phic", phi);
    surfaceScalarField phir("phir", phia - phib);




    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        // Create the limiter to be used for all phase-fractions
        scalarField allLambda(mesh.nFaces(), 1.0);

        // Split the limiter into a surfaceScalarField
        slicedSurfaceScalarField lambda
        (
            IOobject
            (
                "lambda",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimless,
            allLambda,
            false   // Use slices for the couples
        );


        // Create the complete convection flux for alpha
        surfaceScalarField alphaPhi1 //this is the high order solution + the compression (relative velocity) fluxes
        (
            fvc::flux
            (
                phic,
                alpha,
                scheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, beta, schemer),
                alpha,
                schemer
            )
        );

        // Create the bounded (upwind) flux for alpha
        surfaceScalarField alphaPhi1BD //this is the low order solution  without the compression (relative velocity) flux
        (
            upwind<scalar>(mesh, phic).flux(alpha)
        );

        // Calculate the flux correction for alpha
        alphaPhi1 -= alphaPhi1BD;

        // Calculate the limiter for alpha1
        //if (LTS)
        //{
            //const volScalarField& rDeltaT =
                //fv::localEulerDdt::localRDeltaT(mesh);

            //MULES::limiter
            //(
                //allLambda,
                //rDeltaT,
                //geometricOneField(),
                //alpha,
                //alphaPhi1BD,
                //alphaPhi1,
                //zeroField(),
                //zeroField(),
                //oneField(),
                //zeroField()
            //);
        //}
        //else
        //{
            MULES::limiter
            (
                allLambda,
                1.0/runTime.deltaT().value(),
                geometricOneField(),
                alpha,
                alphaPhi1BD,
                alphaPhi1,
                zeroField(),
                zeroField(),
                UniformField<scalar>(0.99*alphaMax.value()),
                zeroField()
            );
        //}

        alphaPhi1 = alphaPhi1BD + lambda*alphaPhi1;

        // Reset allLambda to 1.0
        allLambda = 1.0;

       

        // Solve for alpha1
        solve(fvm::ddt(alpha) + fvc::div(alphaPhi1));


        beta = 1.0 - alpha;

 




        //// Limit to 0.99*alphaMax to avoid extremely small time steps

//#if (defined(OPENFOAM) && (OPENFOAM >= 1906)) || \
    //(defined(OFVERSION) && (OFVERSION >= 600))
    ///* >= OpenFOAM-v1906 or >= OpenFOAM-6 */

            //MULES::explicitSolve
            //(
                //geometricOneField(),
                //alpha,
                //phi,
                //phiAlpha,
                //zeroField(),
                //zeroField(),
                //UniformField<scalar>(0.99*alphaMax.value()),
                //zeroField()
            //);
//#else
            //MULES::explicitSolve
            //(
                //geometricOneField(),
                //alpha,
                //phi,
                //phiAlpha,
                //zeroField(),
                //zeroField(),
                //0.99*alphaMax.value(),
                //0
            //);
//#endif

        ///*if (packingLimiter)
        //{
            //int i=0;
            //while (i<100 && max(alpha).value()>=alphaMax)
            //{
                //#include "packingLimiter.H"
                //i++;
            //}
        //}*/

        //alpha.min(0.99*alphaMax);
        //alpha.max(0);
        //beta = scalar(1.0) - alpha;



   }
}
    
if (debugInfo)
{
	Info<< "Dispersed phase volume fraction = "
	<< alpha.weightedAverage(mesh.V()).value()
	<< "  Min(alpha) = " << gMin(alpha)
	<< "  Max(alpha) = " << gMax(alpha)
	<< endl;
}

rho = alpha*rhoa + beta*rhob;
