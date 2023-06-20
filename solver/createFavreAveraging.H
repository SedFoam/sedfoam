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

Switch favreAveraging
(
    runTime.controlDict().lookupOrDefault("favreAveraging", false)
);

Switch TKEBudgetLES
(
    runTime.controlDict().lookupOrDefault("TKEBudgetLES", false)
);

if (favreAveraging)
{
    Info<<"Favre-averaging are implemented for constant time-step only"<<endl;
}

scalar avCount = 0;

volVectorField UbPrimeF
(
    IOobject
    (
        "UbPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "UbPrimeF",
        dimVelocity,
        vector (0, 0, 0)
    )
);

volVectorField UaPrimeF
(
    IOobject
    (
        "UaPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "UaPrimeF",
        dimVelocity,
        vector (0, 0, 0)
    )
);

volScalarField pMeanF
(
    IOobject
    (
        "pMeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pMeanF",
        dimPressure,
        0
    )
);

volScalarField pPrimeF
(
    IOobject
    (
        "pPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pPrimeF",
        dimPressure,
        0
    )
);

volScalarField RS2
    (
        IOobject
        (
            "RS2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS2",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField RS4
    (
        IOobject
        (
            "RS4",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS4",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField RS2Mean
    (
        IOobject
        (
            "RS2Mean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS2Mean",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField RS4Mean
    (
        IOobject
        (
            "RS4Mean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS4Mean",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volTensorField alphaSqrUaMean
(
    IOobject
    (
        "alphaSqrUaMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "alphaSqrUaMean",
        dimVelocity*dimVelocity,
        tensor::zero
    )
);

volTensorField betaSqrUbMean
(
    IOobject
    (
        "betaSqrUbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "betaSqrUbMean",
        dimVelocity*dimVelocity,
        tensor::zero
    )
);

volVectorField betaUbMean
    (
        IOobject
        (
            "betaUbMean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "betaUbMean",
            dimVelocity,
            vector (0, 0, 0)
        )
    );

volVectorField alphaUaMean
    (
        IOobject
        (
            "alphaUaMean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "alphaUaMean",
            dimVelocity,
            vector (0, 0, 0)
        )
    );

volScalarField alpha_aMean
    (
        IOobject
        (
            "alpha_aMean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "alphaUaMean",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.0
        )
    );

volVectorField UbMeanF
(
    IOobject
    (
        "UbMeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "UbMeanF",
        dimVelocity,
        vector (0, 0, 0)
    )
);

volVectorField UaMeanF
(
    IOobject
    (
        "UaMeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "UaMeanF",
        dimVelocity,
        vector (0, 0, 0)
    )
);
volTensorField UbPrime2MeanF
(
    IOobject
    (
        "UbPrime2MeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "UbPrime2MeanF",
        dimVelocity*dimVelocity,
        tensor::zero
    )
);

volTensorField UaPrime2MeanF
(
    IOobject
    (
        "UaPrime2MeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "UaPrime2MeanF",
        dimVelocity*dimVelocity,
        tensor::zero
    )
);

volScalarField TKEMean
(
    IOobject
    (
        "TKEMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TKEMean",
        dimVelocity*dimVelocity,
        0
    )
);
volScalarField TKEMeanProd
(
    IOobject
    (
        "TKEMeanProd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TKEMeanProd",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0
    )
);

volSymmTensorField fluctstrainTensor
(
    IOobject
    (
        "fluctstrainTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "fluctstrainTensor",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);
volSymmTensorField strainTensor
(
    IOobject
    (
        "strainTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "strainTensor",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);
volScalarField turbDiffusionMean
(
    IOobject
    (
        "turbDiffusionMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "turbDiffusionMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);
volVectorField turbDiffMean
(
    IOobject
    (
        "turbDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "turbDiffMean",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);
volScalarField pressDiffusionMean
(
    IOobject
    (
        "pressDiffusionMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pressDiffusionMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);
volVectorField pressDiffMean
(
    IOobject
    (
        "pressDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "pressDiffMean",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);
volScalarField viscDissMean
(
    IOobject
    (
        "viscDissMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "viscDissMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);
volScalarField viscDiffusionMean
(
    IOobject
    (
        "viscDiffusionMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "viscDiffusionMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);
volScalarField SGSDissMean
(
    IOobject
    (
        "SGSDissMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDissMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);
volVectorField SGSDiffMean
(
    IOobject
    (
        "SGSDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "SGSDiffMean",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);
volScalarField SGSDiffusionMean
(
    IOobject
    (
        "SGSDiffusionMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDiffusionMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);
volScalarField dragContbMean
(
    IOobject
    (
        "dragContbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "dragContbMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volVectorField betaFluct
(
    IOobject
    (
        "betaFluct",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "betaFluct",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volVectorField betaFluctMean
(
    IOobject
    (
        "betaFluctMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "betaFluctMean",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

// Checks if field betaFluctMean is present
IOobject betaFluctMeanHeader
    (
        "betaFluctMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

if (betaFluctMeanHeader.typeHeaderOk<volVectorField>(true))
{
    volVectorField betaFluctMean2
        (
            IOobject
            (
                "betaFluctMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
    betaFluctMean = betaFluctMean2;
}

// Checks if field alpha_aMean is present
IOobject alpha_aMeanHeader
    (
        "alpha_aMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

if (alpha_aMeanHeader.typeHeaderOk<volScalarField>(true))
{
    avCount = runTime.timeIndex();
    Info << avCount << endl;
    volScalarField alpha_aMean2
        (
            IOobject
            (
                "alpha_aMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
    alpha_aMean = alpha_aMean2;
    // Checks if field RS2Mean is present
    IOobject RS2MeanHeader
    (
        "RS2Mean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (RS2MeanHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField RS2Mean2
        (
            IOobject
            (
                "RS2Mean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        RS2Mean = RS2Mean2;
    }
    // Checks if field RS4Mean is present
    IOobject RS4MeanHeader
    (
        "RS4Mean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (RS4MeanHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField RS4Mean2
        (
            IOobject
            (
                "RS4Mean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        RS4Mean = RS4Mean2;
    }
    // Checks if field UbMeanF is present
    IOobject UbMeanFHeader
    (
        "UbMeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (UbMeanFHeader.typeHeaderOk<volVectorField>(true))
    {
        volVectorField UbMeanF2
        (
            IOobject
            (
                "UbMeanF",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        UbMeanF = UbMeanF2;
        betaUbMean = UbMeanF*(1-alpha_aMean);
    }

    // Checks if field UaMeanF is present
    IOobject UaMeanFHeader
    (
        "UaMeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (UaMeanFHeader.typeHeaderOk<volVectorField>(true))
    {
        volVectorField UaMeanF2
        (
            IOobject
            (
                "UaMeanF",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        UaMeanF = UaMeanF2;
        alphaUaMean = UaMeanF*(alpha_aMean+1e-20);
    }

    // Checks if field UbPrime2MeanF is present
    IOobject UbPrime2MeanFHeader
    (
        "UbPrime2MeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (UbPrime2MeanFHeader.typeHeaderOk<volTensorField>(true))
    {
        volTensorField UbPrime2MeanF2
        (
            IOobject
            (
                "UbPrime2MeanF",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        UbPrime2MeanF = UbPrime2MeanF2;
        betaSqrUbMean = (UbPrime2MeanF+UbMeanF*UbMeanF)*(1-alpha_aMean);
    }

    // Checks if field UaPrime2MeanF is present
    IOobject UaPrime2MeanFHeader
    (
        "UaPrime2MeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (UaPrime2MeanFHeader.typeHeaderOk<volTensorField>(true))
    {
        volTensorField UaPrime2MeanF2
        (
            IOobject
            (
                "UaPrime2MeanF",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        UaPrime2MeanF = UaPrime2MeanF2;
        alphaSqrUaMean = (UaPrime2MeanF+UaMeanF*UaMeanF)
                        *(alpha_aMean+1e-20);
    }
}
