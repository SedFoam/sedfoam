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
        alpha.correctBoundaryConditions();
        fvScalarMatrix alphaEqn
        (
             fvm::ddt(alpha)
           + 1/(1-alphaVeg)*fvm::div(phic, alpha, scheme)
           + 1/(1-alphaVeg)*fvm::div(-fvc::flux(-phir, (1.0-alpha-alphaVeg), schemer), alpha, schemer)
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

        beta = scalar(1.0) - alpha - alphaVeg;

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
