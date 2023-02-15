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

#include "phaseModel.H"
#include "fixedValueFvPatchFields.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& phaseName
)
:
    dict_
    (
        transportProperties.subDict("phase" + phaseName)
    ),
    name_(phaseName),
    d_
    (
        dimensionedScalar::getOrDefault
        (
            "d",
            dict_,
            dimensionSet(0, 1, 0, 0, 0, 0, 0),
            200e-6
        )
    ),
    aE_
    (
        dimensionedScalar::getOrDefault
        (
            "aE",
            dict_,
            dimless,
            780
        )
    ),
    bE_
    (
        dimensionedScalar::getOrDefault
        (
            "aE",
            dict_,
            dimless,
            1.8
        )
    ),
    sF_
    (
        dimensionedScalar::getOrDefault
        (
            "sF",
            dict_,
            dimless,
            1
        )
    ),
    hExp_
    (
        dimensionedScalar::getOrDefault
        (
            "hExp",
            dict_,
            dimless,
            2.65
        )
    ),
    nu_
    (
        dimensionedScalar::getOrDefault
        (
            "nu",
            dict_,
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            1e-6
        )
    ),
    rho_
    (
        dimensionedScalar::getOrDefault
        (
            "rho",
            dict_,
            dimensionSet(1, -3, 0, 0, 0, 0, 0),
            1e3
        )
    ),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alpha_
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha", dimless, 0)
    )
{
    const word phiName = IOobject::groupName("phi", phaseName);

    IOobject phiHeader
    (
        phiName,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ
    );

    if (phiHeader.typeHeaderOk<surfaceScalarField>(true))
//    if (phiHeader.headerOk())
    {
        Info<< "Reading face flux field " << phiName << endl;

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U_.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U_.boundaryField(), i)
        {
            if (isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i]))
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(U_) & mesh.Sf(),
                phiTypes
            )
        );
    }
}


Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::New
(
    const fvMesh& mesh,
    const dictionary& transportProperties,
    const word& phaseName
)
{
    return autoPtr<phaseModel>
    (
        new phaseModel(mesh, transportProperties, phaseName)
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// ************************************************************************* //
