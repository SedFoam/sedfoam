/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "sedIncompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::sedIncompressibleTurbulenceModel<TransportModel>::
sedIncompressibleTurbulenceModel
(
    const word& type,
    const volScalarField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transportModel,
    const word& propertiesName
)
:
    TurbulenceModel
    <
        volScalarField,
        geometricOneField,
        incompressibleTurbulenceModel,
        TransportModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transportModel,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::sedIncompressibleTurbulenceModel<TransportModel>>
Foam::sedIncompressibleTurbulenceModel<TransportModel>::New
(
    const volScalarField& alpha,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transportModel,
    const word& propertiesName
)
{
    return autoPtr<sedIncompressibleTurbulenceModel>
    (
        static_cast<sedIncompressibleTurbulenceModel*>(
        TurbulenceModel
        <
            volScalarField,
            geometricOneField,
            incompressibleTurbulenceModel,
            TransportModel
        >::New
        (
            alpha,
            geometricOneField(),
            U,
            alphaRhoPhi,
            phi,
            transportModel,
            propertiesName
        ).ptr())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransportModel>
Foam::tmp<Foam::volScalarField>
Foam::sedIncompressibleTurbulenceModel<TransportModel>::pPrime() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("pPrime", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimPressure, Zero)
    );
}


template<class TransportModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::sedIncompressibleTurbulenceModel<TransportModel>::pPrimef() const
{
    return tmp<surfaceScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("pPrimef", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimPressure, Zero)
    );
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::sedIncompressibleTurbulenceModel<TransportModel>::devReff() const
{
    return devRhoReff();
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::sedIncompressibleTurbulenceModel<TransportModel>::devReff
(
    const volVectorField& U
) const
{
    return devRhoReff(U);
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::sedIncompressibleTurbulenceModel<TransportModel>::divDevReff
(
    volVectorField& U
) const
{
    return divDevRhoReff(U);
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::sedIncompressibleTurbulenceModel<TransportModel>::devRhoReff() const
{
    NotImplemented;

    return devReff();
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::sedIncompressibleTurbulenceModel<TransportModel>::devRhoReff
(
    const volVectorField& U
) const
{
    NotImplemented;

    return nullptr;
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::sedIncompressibleTurbulenceModel<TransportModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    NotImplemented;

    return divDevReff(U);
}

template<class TransportModel>
Foam::tmp<Foam::volScalarField>
Foam::sedIncompressibleTurbulenceModel<TransportModel>::
spherSigmaSGS()
{
    return spherSigmaSGS();
}
// ************************************************************************* //
