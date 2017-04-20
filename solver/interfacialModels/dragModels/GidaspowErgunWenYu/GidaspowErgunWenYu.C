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

#include "GidaspowErgunWenYu.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GidaspowErgunWenYu, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        GidaspowErgunWenYu,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GidaspowErgunWenYu::GidaspowErgunWenYu
(
    const dictionary& interfaceDict,
    const volScalarField& alpha,
    const phaseModel& phasea,
    const phaseModel& phaseb
)
:
    dragModel(interfaceDict, alpha, phasea, phaseb)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GidaspowErgunWenYu::~GidaspowErgunWenYu()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GidaspowErgunWenYu::K
(
    const volScalarField& Ur
) const
{
    volScalarField beta = max(scalar(1) - alpha_, scalar(1.0e-6));

    //    volScalarField bp = pow(beta, -2.65);
    volScalarField bp = pow(beta, -phasea_.hExp());
    volScalarField Re = max(beta*Ur*phasea_.d()*phasea_.sF()/phaseb_.nu(), scalar(1.0e-12));

    volScalarField Cds = 24.0*(1.0 + 0.15*pow(Re, 0.687))/Re;

    forAll(Re, celli)
    {
        if (Re[celli] > 1000.0)
        {
            Cds[celli] = 0.44;
        }
    }
    Cds.correctBoundaryConditions();

    // Wen and Yu (1966)
    // modified in May 17, 2012, by C.Z.
    tmp<volScalarField> tKWenYu = (0.75*Cds*phaseb_.rho()*Ur*bp/(phasea_.d()*phasea_.sF()));
    volScalarField& KWenYu = tKWenYu();
    
    // Ergun
    forAll(beta, cellj)
    {
         if (beta[cellj] < 0.8)
         {
              KWenYu[cellj] =
                      150.0*alpha_[cellj]*phaseb_.nu().value()*phaseb_.rho().value()
                      /sqr(beta[cellj]*phasea_.d().value())
                      + 1.75*phaseb_.rho().value()*Ur[cellj]
                      /(beta[cellj]*phasea_.d().value());
         }
    }
// WARNING: remove this line will makes the parallel computations "instable"
    KWenYu.correctBoundaryConditions();
    
    return tKWenYu;

}

// ************************************************************************* //
