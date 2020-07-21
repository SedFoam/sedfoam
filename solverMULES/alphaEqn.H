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
        #if OFVERSION >= 600
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha,
                phi,
                phiAlpha,
                zeroField(),
                zeroField(),
                UniformField<scalar>(0.99*alphaMax.value()),
                zeroField()
            );
        #elif (OFVERSION < 220 and OFVERSION >= 190) or OFVERSION == 00
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha,
                phi,
                phiAlpha,
                zeroField(),
                zeroField(),
                UniformField<scalar>(0.99*alphaMax.value()),
                zeroField()
            );
        #else
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
        #endif

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
            << "  Min(alpha) = " << gMin(alpha)
            << "  Max(alpha) = " << gMax(alpha)
            << endl;
        }
    }
}

rho = alpha*rhoa + beta*rhob;
