/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "twophaseMixingLength.H"

#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(twophaseMixingLength, 0);
addToRunTimeSelectionTable(RASModel, twophaseMixingLength, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twophaseMixingLength::twophaseMixingLength
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    twophaseRASProperties_
    (
        IOobject
        (
            "twophaseRASProperties",
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    alphaMaxLM_
    (
        twophaseRASProperties_.lookup("alphaMaxLM")
    ),
    kappaLM_
    (
        twophaseRASProperties_.lookup("kappaLM")
    ),
    alpha_(mesh_.lookupObject<volScalarField> ("alpha")),
        k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateEpsilon("epsilon", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
     )
{
  //  nut_ = pow(k_,2)/epsilon_;
  //  nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> twophaseMixingLength::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> twophaseMixingLength::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
	    -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}

tmp<fvVectorMatrix> twophaseMixingLength::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> twophaseMixingLength::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}
bool twophaseMixingLength::read()
{
    if (RASModel::read())
    {
      return true;
    }
    else
    {
      return false;
    }
}
  
void twophaseMixingLength::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

//
//
// Mixing Lenght turbulence model (only for 1D cases with Y the wall-normal direction)
//
//
    dimensionedScalar LmSmall
      (
       "LmSmall",
       dimensionSet(0, 1, 0, 0, 0, 0, 0),
       scalar(1e-4)
       );
    
    //Efective fluid viscosity
    tmp<volTensorField> gradU = fvc::grad(U_);
    volSymmTensorField D = symm(gradU);
    volScalarField magD = ::sqrt(2.0)*mag(D);	
    gradU.clear();

    volVectorField centres = U_.mesh().C();
    volScalarField Y = centres.component(1);
    scalarField& alphaCells = alpha_.internalField();
    scalarField& magDCells = magD.internalField();
    scalarField& YCells = Y.internalField();

    scalar Lm = 0.;
    scalar dY;
    scalar kappaLMs=kappaLM_.value();
    scalar alphaMaxLMs = alphaMaxLM_.value();
    scalar LmPhi = 0.;
    scalar cmu34 = 0.1643;

    nut_.storePrevIter();

    forAll(YCells, cellI)
      {
	dY = YCells[cellI]-YCells[cellI-1];
	LmPhi = LmPhi
	  + kappaLMs*max(scalar(1.0) - Foam::pow(min(alphaCells[cellI]/alphaMaxLMs,scalar(1.0))
						 ,scalar(1.66)), 0.)*dY;
	Lm = LmPhi;
	nut_[cellI] = pow(Lm,2)*magDCells[cellI];
	// for kinetic theory k is required
	k_[cellI] = pow(nut_[cellI]/(0.8*max(Lm,1e-4)),2);
	epsilon_[cellI] = cmu34*Foam::pow(k_[cellI],1.5)/max(Lm,1e-4);
      }
    nut_.relax();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
