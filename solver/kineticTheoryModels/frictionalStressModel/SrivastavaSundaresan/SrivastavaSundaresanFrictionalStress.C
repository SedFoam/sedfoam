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

#include "SrivastavaSundaresanFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SrivastavaSundaresanFrictionalStress, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        SrivastavaSundaresanFrictionalStress,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SrivastavaSundaresanFrictionalStress::SrivastavaSundaresanFrictionalStress
(
    const dictionary& dict
)
:
    frictionalStressModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SrivastavaSundaresanFrictionalStress::~SrivastavaSundaresanFrictionalStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::SrivastavaSundaresanFrictionalStress::
frictionalPressure
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta,
    const dimensionedScalar& p
) const
{

    return
        Fr*pow(max(alpha - alphaMinFriction, scalar(0)), eta)
       /pow(max(alphaMax - alpha, scalar(1.0e-20)), p);
}


Foam::tmp<Foam::volScalarField> Foam::SrivastavaSundaresanFrictionalStress::
frictionalPressurePrime
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta,
    const dimensionedScalar& p
) const
{
    return Fr*
    (
        eta*pow(max(alpha - alphaMinFriction, scalar(0)), eta - 1.0)
       *(alphaMax-alpha) + p*pow(max(alpha - alphaMinFriction, scalar(0)), eta)
    )/pow(max(alphaMax - alpha, scalar(1.0e-20)), p + 1.0);
}

Foam::tmp<Foam::volScalarField> Foam::SrivastavaSundaresanFrictionalStress::muf
(
    const volScalarField& alpha,
    const volScalarField& Theta,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const volScalarField& pf,
    const volSymmTensorField& D,
    const dimensionedScalar& phi
) const
{
    const scalar I2Dsmall = 1.0e-35;

    //Creating muf assuming it should be 0 on the boundary which may not be
    // true
       volScalarField tmuf
       (
          IOobject
          (
             "muf",
             alpha.mesh().time().timeName(),
             alpha.mesh()
          ),
          alpha.mesh(),
          dimensionedScalar("muf", dimensionSet(1, -1, -1, 0, 0), 1e00)
    );
    
    volScalarField& muff = tmuf;
   
      forAll(D, celli)
    {
        if (alpha[celli] >= alphaMinFriction.value())
        {
            muff[celli] =
                0.5*pf[celli]*sin(phi.value())
               /(
                    sqrt(1.0/6.0*(sqr(D[celli].xx() - D[celli].yy())
                  + sqr(D[celli].yy() - D[celli].zz())
                  + sqr(D[celli].zz() - D[celli].xx()))
                  + sqr(D[celli].xy()) + sqr(D[celli].xz())
                  + sqr(D[celli].yz())) + I2Dsmall
                );
        }
        if (alpha[celli] < alphaMinFriction.value())
        {
            muff[celli] = 0.0;
        }
    }

    muff.correctBoundaryConditions();

    return tmuf;
}
 
// ************************************************************************* //
