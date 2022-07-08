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
if (this->runTime_.outputTime())
{
    volScalarField Kprod
    (
       IOobject
       (
            "Kprod",
            k_.time().timeName(),
            k_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       G
    );
    volScalarField Kdiff
    (
       IOobject
       (
            "Kdiff",
            k_.time().timeName(),
            k_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       fvc::laplacian(DkEff(), k_)
    );
    volScalarField Kdissip
    (
       IOobject
       (
            "Kdissip",
            k_.time().timeName(),
            k_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       -epsilon_
    );
    
    volScalarField Kdissdrag1
    (
       IOobject
       (
            "Kdissdrag1",
            k_.time().timeName(),
            k_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       ESD_*k_
    );
    volScalarField Kdissdrag2
    (
       IOobject
       (
            "Kdissdrag2",
            k_.time().timeName(),
            k_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       KE4_*ESD4_*nut
    );
    volScalarField Kdissdrag3
    (
       IOobject
       (
            "Kdissdrag3",
            k_.time().timeName(),
            k_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       ESD2()*KE2_ * k_
    );

    this->runTime_.write();

}
