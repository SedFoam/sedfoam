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
    word schemeRey("laplacian(nut,alpha)");

    surfaceScalarField phic("phic", phi);
    surfaceScalarField phir("phir", phia - phib);




    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        alpha.correctBoundaryConditions();

        fvScalarMatrix alphaEqn
        (
             fvm::ddt(alpha)
           + fvm::div(phic, alpha, scheme)
           + fvm::div(-fvc::flux(-phir, (1.0-alpha), schemer), alpha, schemer)
           ==
           //fvm::laplacian(turbulenceb->nut(),alpha)
           //pos(alphaPlastic-alpha)*fvm::laplacian(turbulenceb->nut(),alpha)
           fvc::div(fvc::interpolate(turbulenceb->nut())*mesh.magSf()*fvc::snGrad(alpha))
           
        );    
        alphaEqn.relax();
        alphaEqn.solve();

        if (packingLimiter)
        {
            int i=0;
            while (i<100 && max(alpha).value()>=alphaMax.value())
            {
                #include "packingLimiter.H"
                i++;
            }
        }
        alpha.min(alphaMax);
        alpha.max(0);

        beta = scalar(1.0) - alpha;

        if (debugInfo)
        {
            Info<< "Dispersed phase volume fraction = "
            << alpha.weightedAverage(mesh.V()).value()
            << "  Min(alpha) = " << gMin(alpha)
            << "  Max(alpha) = " << gMax(alpha)
            << endl;
        }
    }
}

rho = alpha*rhoa + beta*rhob;

if (bidispersity)
{
	
	
		//volScalarField ALPHA_s=alpha_s/(alpha_s+alpha_l+alphaSmall);
		//volScalarField ALPHA_l=alpha_l/(alpha_s+alpha_l+alphaSmall);
		volScalarField d_m=ALPHA_l*d_l+ALPHA_s*d_s;

		gradUaT = fvc::grad(Ua)().T();
		volScalarField magD_sl
		(
			::sqrt(2.0)*mag(symm(gradUaT))
		);


		dimensionedScalar Beta_sl
		(
		"Beta_sl",
		dimensionSet(0, 0, 0, 0, 0, 0, 0),
		0.3744
		);

		dimensionedScalar C_sl
		(
		"C_sl",
		dimensionSet(0, 0, 0, 0, 0, 0, 0),
		0.2712
		);

		dimensionedScalar E_sl
		(
		"E_sl",
		dimensionSet(0, 0, 0, 0, 0, 0, 0),
		2.0957
		);

		dimensionedScalar A_sl
		(
		"A_sl",
		dimensionSet(0, 0, 0, 0, 0, 0, 0),
		0.108
		);
		dimensionedScalar dsmall
		(
		"dsmall",
		dimensionSet(0, 1, 0, 0, 0, 0, 0),
		1e-8
		);		
		word scheme("div(phi,alpha)");
	    surfaceScalarField phic("phic", phi);

		dimensionedScalar R_sl=d_l/d_s;
		volScalarField Lambda_sl=1-ALPHA_s;
		F=(R_sl-1)+E_sl*Lambda_sl*sqr(R_sl-1);
		volScalarField F=(R_sl-1)+E_sl*Lambda_sl*sqr(R_sl-1);
		//volScalarField F=Lambda_sl*scalar(0)+1;
	

    dimensionedScalar dproves
    (
        "dproves",
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        0.1
    );
    
    

    dimensionedScalar MagDsmall
    (
        "MagDsmall",
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        1e-19
    );
    
    dimensionedScalar alphaCheck
    (
        "alphaCheck",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        0.00001
    );
    dimensionedVector gravityVertical
    (
        "gravityVertical",
        dimensionSet(0, 1, -2, 0, 0, 0, 0),
        vector(0,-8.89,0)
    );
		//f_sl=pos(alpha-alphaCheck)*ALPHA_l*(Beta_sl*rhoa*mag(gravity)*pos(magD_sl-MagDsmall)*dproves*sqr(d_m)*F/(C_sl*rhoa*mag(gravity)*max(d_m,dsmall)+pS))*gravity/mag(gravity);
		f_sl=pos(alpha-alphaCheck)*ALPHA_l*(Beta_sl*rhoa*mag(gravity)*pos(magD_sl-MagDsmall)*dproves*sqr(d_m)*F/(C_sl*rhoa*mag(gravity)*max(d_m,dsmall)+0*pS))*gravity/mag(gravity);

	    phiG =  linearInterpolate(f_sl) & mesh.Sf();
		//surfaceScalarField phiG =  linearInterpolate(f_sl) & mesh.Sf();

		volScalarField D_sl=A_sl*magD_sl*sqr(d_m);

        fvScalarMatrix alphaBiEqn
        (
             fvm::ddt(ALPHA_s)
           + fvm::div(phic, ALPHA_s, scheme)
           + fvm::div(phiG,ALPHA_s) //segregation term
         //  + fvc::div(phiG) //segregation term
         //  - fvm::laplacian(D_sl, ALPHA_s)   // diffusion term 
           
        );

        alphaBiEqn.relax();
        alphaBiEqn.solve();
        
        ALPHA_s.min(1);
        ALPHA_s.max(0);
        
        ALPHA_s=pos(alpha-alphaCheck)*ALPHA_s;
        ALPHA_l=pos(alpha-alphaCheck)*(1-ALPHA_s);
        
        if (debugInfo)
        {
            Info<< "Large particle phase volume fraction = "
            << ALPHA_l.weightedAverage(mesh.V()).value()
            << "  Min(ALPHA_l) = " << gMin(ALPHA_l)
            << "  Max(ALPHA_l) = " << gMax(ALPHA_l)
            << endl;
            
            Info<< "Small particle phase volume fraction = "
            << ALPHA_s.weightedAverage(mesh.V()).value()
            << "  Min(ALPHA_s) = " << gMin(ALPHA_s)
            << "  Max(ALPHA_s) = " << gMax(ALPHA_s)
            << endl;
            
        }

}	
