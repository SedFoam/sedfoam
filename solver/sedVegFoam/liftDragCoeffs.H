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
volScalarField magUr
(
    mag(-Ua + Ub)
);

volScalarField magUb
(
    mag(Ub)
);

volScalarField Ka(draga->K(magUr));

K = Ka;

if (dragPhase == "b")
{
    volScalarField Kb(dragb->K(magUr));
    K = Kb;
}
else if (dragPhase == "blended")
{
    volScalarField Kb(dragb->K(magUr));
    K = (beta*Ka + alpha*Kb);
}

Kveg = dragbveg->K(magUb);

volVectorField liftCoeff
(
    Cl*(beta*rhob + alpha*rhoa)*((-Ua + Ub) ^ fvc::curl(U))
);
