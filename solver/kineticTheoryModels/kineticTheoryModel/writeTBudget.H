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
if (Ua_.time().outputTime())
{
    volScalarField psTerm
    (
       IOobject
       (
            "psTerm",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       -fvc::SuSp((PsCoeff*I) && dU, Theta_)
    );

    volScalarField Tprod
    (
       IOobject
       (
            "Tprod",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       tau && dU
    );
    volScalarField Tdiff
    (
       IOobject
       (
            "Tdiff",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       fvc::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
    );
    volScalarField Tdissip
    (
       IOobject
       (
            "Tdissip",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       fvc::Sp(-gammaCoeff, Theta_)
    );
    volScalarField TdragJ1
    (
       IOobject
       (
            "TdragJ1",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       fvc::Sp(-J1, Theta_)
    );
    volScalarField TdragJ2
    (
       IOobject
       (
            "TdragJ2",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       (1-killJ2_)*J2
    );
    Ua_.time().write();
}
