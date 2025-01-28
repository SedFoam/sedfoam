/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "GidaspowSchillerNaumann_veg.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GidaspowSchillerNaumann_veg, 0);

    addToRunTimeSelectionTable
    (
        vegDragModel,
        GidaspowSchillerNaumann_veg,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GidaspowSchillerNaumann_veg::GidaspowSchillerNaumann_veg
(
    const dictionary& interfaceDict,
    const phaseModel& phasea,
    const phaseModel& phaseb
)
:
    vegDragModel(interfaceDict, phasea, phaseb)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GidaspowSchillerNaumann_veg::~GidaspowSchillerNaumann_veg()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GidaspowSchillerNaumann_veg::K
(
    const volScalarField& Ur
) const
{
    volScalarField beta(max(scalar(1) - phaseb_.alpha(), scalar(1e-6)));
    

    volScalarField bp(pow(beta, -phaseb_.hExp()));

   
    const dimensionedScalar dragConst
    (
        interfaceDict_.getOrDefault
        (
            "dragConst",
            dimensionedScalar("dragConst",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          1)
        )
    );

    const dimensionedScalar Coe = dragConst;
   
	   Info<<"dragConst in GidaspowSchiller = \t"<<dragConst<<endl; 

    volScalarField Re
    (
        max(beta*Ur*phaseb_.d()*phaseb_.sF()/phasea_.nu(), scalar(1.0e-9))
    );

    volScalarField Cds
    (
        neg(Re - 1000)*(24.0*(1.0 + 0.15*pow(Re, 0.687))/Re)
      + pos(Re - 1000)*0.44
    );

    return 0.75*Coe*Cds*phasea_.rho()*Ur*bp/(phaseb_.d()*phaseb_.sF());
}


// ************************************************************************* //
