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

#include "SchillerNaumann_veg.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SchillerNaumann_veg, 0);

    addToRunTimeSelectionTable
    (
        vegDragModel,
        SchillerNaumann_veg,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SchillerNaumann_veg::SchillerNaumann_veg
(
    const dictionary& interfaceDict,
    const phaseModel& phasea,
    const phaseModel& phaseb
)
:
    vegDragModel(interfaceDict, phasea, phaseb)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SchillerNaumann_veg::~SchillerNaumann_veg()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::SchillerNaumann_veg::K
(
    const volScalarField& Ur
) const
{
    //// Here phase b represents the particle phase //

    volScalarField Re(max(Ur*phaseb_.d()*phaseb_.sF()/
    phaseb_.nu(), scalar(1.0e-3)));

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


    volScalarField Cds
    (
        neg(Re - 1000)*(24.0*(1.0 + 0.15*pow(Re, 0.687))/Re)
      + pos(Re - 1000)*0.44
    );

    return 0.75*Coe*Cds*phasea_.rho()*Ur/(phaseb_.d()*phaseb_.sF());
}

// ************************************************************************* //
