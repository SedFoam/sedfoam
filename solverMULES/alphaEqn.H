{
    word scheme("div(phi,alpha)");
    word schemer("div(phir,alpha)");

    surfaceScalarField phic("phic", phi);
    surfaceScalarField phir("phir", phia - phib);
    
    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        surfaceScalarField phiAlpha
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

        // Limit to 0.99*alphaMax to avoid extremely small time steps
        MULES::explicitSolve
        (
            geometricOneField(),
            alpha,
            phi,
            phiAlpha,
            zeroField(),
            zeroField(),
            0.99*alphaMax,
            0
        );

        /*if (packingLimiter)
        {
	    int i=0;
            while (i<100 && max(alpha).value()>=alphaMax)
            {
                #include "packingLimiter.H"
                i++;
            }
        }*/
        
        alpha.min(0.99*alphaMax);
        alpha.max(0);

        beta = scalar(1.0) - alpha;

	if (debugInfo)
	{
	    Info<< "Dispersed phase volume fraction = "
		<< alpha.weightedAverage(mesh.V()).value()
		<< "  Min(alpha) = " << min(alpha).value()
		<< "  Max(alpha) = " << max(alpha).value()
		<< endl;
        }    
    }
}

rho = alpha*rhoa + beta*rhob;
