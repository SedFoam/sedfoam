/*---------------------------------------------------------------------------*\
Copyright (C) 2015 Cyrille Bonamy, Julien Chauchat, Tian-Jian Hsu
                   and contributors

License
    This file is part of SedFOAM.

    SedFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SedFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with SedFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("eps", dimensionSet(0, 1, -1, 0, 0), 1e-6)
);
volScalarField err
(
    IOobject
    (
        "err",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("err", dimensionSet(0, 1, -1, 0, 0), 0.1)
);

volScalarField Ufalltmp
(
    IOobject
    (
        "Ufalltmp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("small_", dimensionSet(0, 1, -1, 0, 0), 1.e-6)
);


volScalarField UfallField
(
    IOobject
    (
        "UfallField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimensionSet(0, 1, -1, 0, 0), 1.0)
);

scalar iter=0;
while ((err > eps) and (iter < 10))
{
    iter = iter + 1;
    UfallField = mag(gravity)*(rhoa-rhob)/draga->K(Ufalltmp);
    err = max(mag(UfallField - Ufalltmp));
    // under-relaxation for convergence stability
    Ufalltmp = 0.5*Ufalltmp+0.5*UfallField;
    //Info<<iter<<endl;
}

if (iter == 10)
{
    Info << "Maximum number of iterations reached for settling "
         << "velocity calculation"<<endl;
}
Info << "Settling velocity in quiescent water: Ufall="
     << max(UfallField.internalField()).value()
     << " with residual:" << max(err).value() << endl;
