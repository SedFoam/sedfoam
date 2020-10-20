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
volTensorField gradUaT = fvc::grad(Ua)().T();
volTensorField gradUbT = fvc::grad(Ub)().T();
volSymmTensorField devS = dev(symm(gradUaT));

if (kineticTheory.on())
{
    // Compute Kinetic Theory including granular temperature solution
    kineticTheory.solve
    (
        gradUaT, turbulenceb->k(), turbulenceb->epsilon(),
        turbulenceb->nut(), B, runTime
    );

    // Update solid phase viscosities and collisional pressure
    nuEffa = kineticTheory.mua()/((alpha + alphaSmall)*rhoa);
    lambdaUa = kineticTheory.lambda();
    pa = kineticTheory.pa();

    if (debugInfo)
    {
        Info<< "Contact pressure  pff: Min =" << gMin(pff)
            <<", Max =" << gMax(pff)<<endl;
        Info<< "Collisional press. pa: Min =" << gMin(pa)
            <<", Max =" << gMax(pa)<<endl;
    }
}
else               // If not using kinetic theory is using Ct model
{
    nuEffa = sqr(Ct)*turbulenceb->nut() + nua;
}
