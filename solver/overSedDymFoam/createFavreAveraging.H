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
//##include "getTimeIndex.H"
// #include <cstdlib>

/* Variables defined for Fluid part in the TKE budget for intiation  */

dimensionedScalar zeroo
(
    "zeroo",
    dimensionSet(0, 2, -2, 0, 0, 0, 0),
    scalar(0.0)
);
Switch favreAveraging
(
    runTime.controlDict().lookupOrDefault("favreAveraging", false)
);

Switch favreAverage_fluid
(
    runTime.controlDict().lookupOrDefault("favreAverage_fluid", false)
);

Switch favreAverage_solid
(
    runTime.controlDict().lookupOrDefault("favreAverage_solid", false)
);

Switch TKEbudgetLES_fluid
(
    runTime.controlDict().lookupOrDefault("TKEbudgetLES_fluid", false)
);

if (TKEbudgetLES_fluid == true)
{
    favreAverage_fluid = true;
    Info<<"turn favreAverage_fluid to true for TKEbudgetLES_fluid"<<endl;
}

Switch TKEbudgetLES_solid
(
    runTime.controlDict().lookupOrDefault("TKEbudgetLES_solid", false)
);

if (TKEbudgetLES_solid == true)
{
    favreAverage_solid = true;
    Info<<"turn favreAverage_solid to true for TKEbudgetLES_solid"<<endl;
}

Switch Reset_Average
(
  runTime.controlDict().lookupOrDefault("Reset_Average", false)
);

if (Reset_Average)
{
  Info<<"Reset fluid options = \t"<<Reset_Average<<endl;
}

Info<<"Reset fluid Mode = \t"<<Reset_Average<<endl;

scalar time_average = 0;
scalar deltaT = 0;
scalar T0 = 0;

if (favreAveraging)
{
  T0 =runTime.controlDict().get<scalar>("StartAverageTime");
}

Info<<"Start average time: T0 = \t"<<T0<<endl;

volScalarField dragVegMean
(
    IOobject
    (
        "dragVegMean",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "dragVegMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField WbcubeMean
(
    IOobject
    (
        "WbcubeMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "WbcubeMean",
        dimVelocity*dimVelocity*dimVelocity,
        0.0
    )
);

volScalarField UbSqrWbMean
(
    IOobject
    (
        "UbSqrWbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "UbSqrWbMean",
        dimVelocity*dimVelocity*dimVelocity,
        0.0
    )
);

volVectorField turbDiffusionMeanWallNormal_b
(
    IOobject
  (
    "turbDiffusionMeanWallNormal_b",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
  mesh,
  dimensionedVector
    (
        "turbDiffusionMeanWallNormal_b",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        vector::zero
    )
);

volScalarField TKEMean_b
(
    IOobject
    (
        "TKEMean_b",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TKEMean_b",
        dimVelocity*dimVelocity,
        0
    )
);

volScalarField TKEMeanProd_b
(
    IOobject
    (
        "TKEMeanProd_b",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TKEMeanProd_b",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0
    )
);

volSymmTensorField fluctstrainTensor_b
(
    IOobject
    (
        "fluctstrainTensor_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "fluctstrainTensor_b",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);

volSymmTensorField strainTensor_b
(
    IOobject
    (
        "strainTensor_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "strainTensor_b",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);

volScalarField viscDissMeanI_b
(
    IOobject
    (
        "viscDissMeanI_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "viscDissMeanI_b",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField RS2_b
    (
        IOobject
        (
            "RS2_b",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS2_b",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField RS4_b
    (
        IOobject
        (
            "RS4_b",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS4_b",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField pMean
(
    IOobject
    (
        "pMean",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pMean",
        dimPressure,
        0
    )
);

volTensorField UbPrime2MeanF
  (
      IOobject
        (
        "UbPrime2MeanF",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
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

volTensorField betaSqrUbMean
(
    IOobject
    (
        "betaSqrUbMean",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
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

volScalarField RS2Mean_b
    (
        IOobject
        (
            "RS2Mean_b",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS2Mean_b",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField RS4Mean_b
    (
        IOobject
        (
            "RS4Mean_b",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS4Mean_b",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField SijbSijb
  (
    IOobject
      (
        "SijbSijb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
    mesh,
    dimensionedScalar
      (
        "SijbSijb",
        dimVelocity/dimLength*dimVelocity/dimLength,
        0
       )
  );

volScalarField SGSDissMean_b
(
    IOobject
    (
        "SGSDissMean_b",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDissMean_b",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField SGSDissMeanI_b
(
    IOobject
    (
        "SGSDissMeanI_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDissMeanI_b",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volSymmTensorField strainTensorM_b
(
    IOobject
    (
        "strainTensorM_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "strainTensorM_b",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);

volScalarField turbDiffusionMean_b
(
    IOobject
    (
        "turbDiffusionMean_b",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "turbDiffusionMean_b",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField pressDiffusionMean_b
(
    IOobject
    (
        "pressDiffusionMean_b",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pressDiffusionMean_b",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField SGSDiffusionMean_b
(
    IOobject
    (
        "SGSDiffusionMean_b",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDiffusionMean_b",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField SijSij_b
  (
    IOobject
      (
        "SijSij_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
    mesh,
    dimensionedScalar
      (
        "SijSij_b",
        dimVelocity/dimLength*dimVelocity/dimLength,
        0
       )
  );

volScalarField viscDiffusionMean_b
(
    IOobject
    (
        "viscDiffusionMean_b",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "viscDiffusionMean_b",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volVectorField SGSDiffMean_b
(
    IOobject
    (
        "SGSDiffMean_b",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "SGSDiffMean_b",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volScalarField dragConstbMean
(
    IOobject
    (
        "dragConstbMean",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "dragConstbMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

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

volScalarField alpha_aMean
    (
    IOobject
        (
      "alpha_aMean",
      runTime.timeName(),
      mesh,
      IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
    mesh,
    dimensionedScalar
      (
      "alpha_aMean",
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
        IOobject::READ_IF_PRESENT,
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

volVectorField betaUbMean
    (
    IOobject
      (
        "betaUbMean",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
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

volScalarField viscDissMean_b
  (
      IOobject
        (
        "viscDissMean_b",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
        ),
      mesh,
      dimensionedScalar
        (
        "viscDissMean_b",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
        )
  );

volScalarField pPrime
(
    IOobject
    (
        "pPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pPrime",
        dimPressure,
        0
    )
);
volScalarField TKEMean_a
(
    IOobject
    (
        "TKEMean_a",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TKEMean_a",
        dimVelocity*dimVelocity,
        0
    )
);

volScalarField TKEMeanProd_a
(
    IOobject
    (
        "TKEMeanProd_a",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TKEMeanProd_a",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0
    )
);


volSymmTensorField fluctstrainTensor_a
(
    IOobject
    (
        "fluctstrainTensor_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "fluctstrainTensor_a",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);

volSymmTensorField strainTensor_a
(
    IOobject
    (
        "strainTensor_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "strainTensor_a",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);

volScalarField viscDissMeanI_a
(
    IOobject
    (
        "viscDissMeanI_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "viscDissMeanI_a",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField RS2_a
    (
        IOobject
        (
            "RS2_a",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS2_a",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField RS4_a
    (
        IOobject
        (
            "RS4_a",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS4_a",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField paMean
(
    IOobject
    (
        "paMean",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "paMean",
        dimPressure,
        0
    )
);

volTensorField UaPrime2MeanF
    (
            IOobject
                (
                    "UaPrime2MeanF",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
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

volTensorField alphaSqrUaMean
(
    IOobject
    (
        "alphaSqrUaMean",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
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

volScalarField RS2Mean_a
    (
        IOobject
        (
            "RS2Mean_a",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS2Mean_a",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField RS4Mean_a
    (
        IOobject
        (
            "RS4Mean_a",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS4Mean_a",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField SijaSija
    (
            IOobject
                (
                    "SijaSija",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
            mesh,
            dimensionedScalar
                (
                    "SijaSija",
                    dimVelocity/dimLength*dimVelocity/dimLength,
                    0
             )
    );

volScalarField SGSDissMean_a
(
    IOobject
    (
        "SGSDissMean_a",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDissMean_a",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField SGSDissMeanI_a
(
    IOobject
    (
        "SGSDissMeanI_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDissMeanI_a",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volSymmTensorField strainTensorM_a
(
    IOobject
    (
        "strainTensor_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "strainTensor_a",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);

volScalarField turbDiffusionMean_a
(
    IOobject
    (
        "turbDiffusionMean_a",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "turbDiffusionMean_a",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField pressDiffusionMean_a
(
    IOobject
    (
        "pressDiffusionMean_a",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pressDiffusionMean_a",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField SGSDiffusionMean_a
(
    IOobject
    (
        "SGSDiffusionMean_a",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDiffusionMean_a",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volVectorField turbDiffMean_a
(
    IOobject
    (
        "turbDiffMean_a",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "turbDiffMean_a",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volScalarField SijSij_a
    (
            IOobject
                (
                    "SijSij_a",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
            mesh,
            dimensionedScalar
                (
                    "SijSij_a",
                    dimVelocity/dimLength*dimVelocity/dimLength,
                    0
             )
    );

volScalarField viscDiffusionMean_a
(
    IOobject
    (
        "viscDiffusionMean_a",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "viscDiffusionMean_a",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
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

volVectorField UaMeanF
    (
            IOobject
                (
                    "UaMeanF",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
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

volVectorField alphaUaMean
        (
            IOobject
                (
                    "alphaUaMean",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
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

volScalarField viscDissMean_a
    (
            IOobject
                (
                    "viscDissMean_a",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
            mesh,
            dimensionedScalar
                (
                    "viscDissMean_a",
                    dimVelocity*dimVelocity*dimVelocity/dimLength,
                    0.0
                )
    );

volScalarField paPrime
(
    IOobject
    (
        "paPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "paPrime",
        dimPressure,
        0
    )
);

/*   Initialization of variables for momentum budget   */

volScalarField pMech_b
(
    IOobject
    (
        "pMech_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pMech_b",
        dimPressure,
        0
    )
);

volScalarField pA
(
    IOobject
    (
        "pA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pA",
        dimPressure,
        0
    )
);

volScalarField pFF
(
    IOobject
    (
        "pFF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pFF",
        dimPressure,
        0
    )
);

volVectorField rho_mixd
(
    IOobject
    (
        "rho_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "rho_mixd",
        dimensionSet(1, -2, -1, 0, 0, 0, 0),
        vector (0, 0, 0)
    )
);

volTensorField nonl_mixd
(
    IOobject
    (
        "nonl_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "nonl_mixd",
        dimensionSet(1, -1, -2, 0, 0, 0, 0),
        tensor::zero
    )
);

volVectorField buoy_mixd
    (
        IOobject
        (
            "buoy_mixd",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "buoy_mixd",
            dimensionSet(1, -2, -2, 0, 0, 0, 0),
            vector (0, 0, 0)
        )
    );

volVectorField SolidPressureContr_b
(
    IOobject
    (
        "SolidPressureContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(p_rbgh)
);

volVectorField divTauaContr_b
(
    IOobject
    (
        "divTauaContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField ViscStressContr_b
(
    IOobject
    (
        "ViscStressContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField FrictContr_b
(
    IOobject
    (
        "FrictContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField phiRaContr_b
(
    IOobject
    (
        "phiRaContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField phiRa2Contr_b
(
    IOobject
    (
        "phiRa2Contr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField divRcaContr_b
(
    IOobject
    (
        "divRcaContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField RcaContr_b
(
    IOobject
    (
        "RcaContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField BulkContr_b
(
    IOobject
    (
        "BulkContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField gravityContr_b
(
    IOobject
    (
        "gravityContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField DragContr_b
(
    IOobject
    (
        "DragContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField SuspContr_b
(
    IOobject
    (
        "SuspContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField ExtContr_b
(
    IOobject
    (
        "ExtContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volScalarField Erosion
        (
            IOobject
                (
                    "Erosion",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
            IOobject::NO_WRITE
                ),
            mesh,
            dimensionedScalar
                (
                    "Erosion",
                        dimensionSet(0, 1, -1, 0, 0, 0, 0),
                        0.0
                )
 );

volScalarField Deposition
        (
            IOobject
                (
                    "Deposition",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
            IOobject::NO_WRITE
                ),
            mesh,
            dimensionedScalar
                (
                    "Deposition",
                        dimensionSet(0, 1, -1, 0, 0, 0, 0),
                        0.0
                )
        );

if (Reset_Average)
{
    time_average = 0;
    Info<<"Reset averaging: time_average= \t"<<time_average<<endl;
}
else
{
    if (runTime.value() < T0)
    {
        Info<<"Averaging will start at t=T0: \t"<<T0<<endl;
    }
    else
    {
        time_average = runTime.value() - T0;
        Info<<"Averaging since: T0= \t"<<T0<<endl;
        Info<<"--- continue averaging at : t   = \t"<<runTime.value()<<endl;
        Info<<"--- with averaging period : t-T0= \t"<<runTime.value()-T0<<endl;
        
    }
}
