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
 #include <cstdlib>

/* Variables defined for Fluid part in the TKE budget for intiation  */

dimensionedScalar zeroo
(
    "zeroo",
    dimensionSet(0, 2, -2, 0, 0, 0, 0),
    scalar(0.0)
);

Switch favreAverage_fluid
(
    runTime.controlDict().lookupOrDefault("favreAverage_fluid", false)
);

Switch AverageMode
(
    runTime.controlDict().lookupOrDefault("AverageMode", false)
);

Switch Reset_fluid
(
  runTime.controlDict().lookupOrDefault("Reset_fluid", false)
);

if (Reset_fluid)
{
  Info<<"Reset fluid options = \t"<<Reset_fluid<<endl;
}

Info<<"Reset fluid Mode = \t"<<Reset_fluid<<endl;

scalar timeaverage_fluid = 0;
double runtime_fluid     = 0;
scalar deltaT = 0;
scalar T0 = 0;
double time_integrate_fluid = 0;
double time_diff_fluid = 0;

if (AverageMode)
{
  T0 =runTime.controlDict().get<scalar>("StartAverageTime");
}

Info<<"To before else looop in Average fluid non reset mode = \t"<<T0<<endl;

volScalarField TKEMeanProdI_b
(
    IOobject
    (
        "TKEMeanProdI_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TKEMeanProdI_b",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0
    )
);

volScalarField dragVegMean
(
    IOobject
    (
        "dragVegMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
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

volVectorField betaFluctMean_b
(
    IOobject
    (
        "betaFluctMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "betaFluctMean_b",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
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

volScalarField UbsquarWbMean
(
    IOobject
    (
        "UbsquarWbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "UbsquarWbMean",
        dimVelocity*dimVelocity*dimVelocity,
        0.0
    )
);

volVectorField turbDiffusionMeanDavid
(
    IOobject
  (
    "turbDiffusionMeanDavid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
  mesh,
  dimensionedVector
    (
        "turbDiffusionMeanDavid",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        vector::zero
    )
);

volScalarField testVegMean_b
  (
  IOobject
    (
    "testVegMean_b",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
        mesh,
        dimensionedScalar
        (
    "testVegMean_b",
    dimensionSet(1, -3, -1, 0, 0, 0, 0),
    0.0
        )
    );

volVectorField UbMean
  (
  IOobject
      (
        "UbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
    mesh,
    dimensionedVector
      (
        "UbMean",
        dimVelocity,
        vector (0, 0, 0)
      )
    );

volScalarField beta_bMean
    (
        IOobject
        (
            "beta_bMean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "beta_bMean",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.0
        )
    );

///
volScalarField TKEMean_b
(
    IOobject
    (
        "TKEMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
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
        IOobject::NO_READ,
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

volVectorField betaFluct_b
(
    IOobject
    (
        "betaFluct_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "betaFluct_b",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
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

volScalarField pMean_b
(
    IOobject
    (
        "pMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pMean_b",
        dimPressure,
        0
    )
);

volTensorField UbPrime2MeanF_b
  (
      IOobject
        (
        "UbPrime2MeanF_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
      mesh,
      dimensionedTensor
        (
        "UbPrime2Mean",
        dimVelocity*dimVelocity,
        tensor::zero
         )
  );

volTensorField betaSqrUbMean_b
(
    IOobject
    (
        "betaSqrUbMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "betaSqrUbMean_b",
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
            IOobject::NO_READ,
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
            IOobject::NO_READ,
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
        IOobject::NO_READ,
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
        IOobject::NO_READ,
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
        IOobject::NO_READ,
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
        IOobject::NO_READ,
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

volVectorField turbDiffMeanJ_b
(
    IOobject
    (
        "turbDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "turbDiffMeanJ_b",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
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
        IOobject::NO_READ,
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
        IOobject::NO_READ,
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

volVectorField viscDiffMeanJ_b
(
    IOobject
    (
        "viscDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "viscDiffMeanJ_b",
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
        IOobject::NO_READ,
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

volVectorField UbPrimeF_b
  (
    IOobject
      (
        "UbPrimeF_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
    mesh,
    dimensionedVector
      (
        "UbPrimeF_b",
        dimVelocity,
        vector (0, 0, 0)
      )
  );

volScalarField alpha_aMean_F
    (
    IOobject
        (
      "alpha_aMean_F",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
          IOobject::NO_WRITE
      ),
    mesh,
    dimensionedScalar
      (
      "alpha_aMean_F",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        0.0
      )
  );

volVectorField UbMeanF_b
  (
      IOobject
        (
        "UbMeanF_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
      mesh,
      dimensionedVector
        (
        "UbMeanF_b",
        dimVelocity,
        vector (0, 0, 0)
        )
    );

volVectorField betaUbMean_F
    (
    IOobject
      (
        "betaUbMean_F",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedVector
        (
          "betaUbMean_F",
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
        IOobject::NO_READ,
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

volScalarField pPrimeF_b
(
    IOobject
    (
        "pPrimeF_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pPrimeF_b",
        dimPressure,
        0
    )
);

/*     Variable for Mass budget terms for initiations */

volVectorField Mean_alphaUa
    (
      IOobject
            (
              "Mean_alphaUa",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
            ),
          mesh,
          dimensionedVector
            (
              "Mean_alphaUa",
              dimVelocity,
              vector (0, 0, 0)
            )
      );

volVectorField Mean_UaFalphaF
    (
      IOobject
            (
              "Mean_UaFalphaF",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
            ),
          mesh,
          dimensionedVector
            (
              "Mean_UaFalphaF",
              dimVelocity,
              vector (0, 0, 0)
            )
      );

volVectorField MeanUaF_MeanalphaF
    (
      IOobject
            (
              "MeanUaF_MeanalphaF",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
            ),
          mesh,
          dimensionedVector
            (
              "MeanUaF_MeanalphaF",
              dimVelocity,
              vector (0, 0, 0)
            )
      );

volVectorField MeanUaF
    (
      IOobject
            (
              "MeanUaF",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
            ),
          mesh,
          dimensionedVector
            (
              "MeanUaF",
              dimVelocity,
              vector (0, 0, 0)
            )
      );

volScalarField alphaFMean
    (
        IOobject
    (
          "alphaFMean",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
      IOobject::NO_WRITE
          ),
        mesh,
        dimensionedScalar
        (
          "alphaFMean",
          dimensionSet(0, 0, 0, 0, 0, 0, 0),
          0.0
        )
    );

volScalarField alphaF
    (
        IOobject
    (
          "alphaF",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
      IOobject::NO_WRITE
          ),
        mesh,
        dimensionedScalar
        (
          "alphaF",
          dimensionSet(0, 0, 0, 0, 0, 0, 0),
          0.0
        )
    );

/*  Variable initialization for particle phase */

volScalarField WacubeMean
(
    IOobject
    (
        "WacubeMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "WacubeMean",
        dimVelocity*dimVelocity*dimVelocity,
        0.0
    )
);

volScalarField alpha_aMean_a
        (
            IOobject
                (
                    "alpha_aMean_a",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
            mesh,
            dimensionedScalar
                (
                    "alpha_aMean_a",
                        dimensionSet(0, 0, 0, 0, 0, 0, 0),
                        0.0
                )
        );

volScalarField UasquarWbMean
(
    IOobject
    (
        "UasquarWbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "UasquarWbMean",
        dimVelocity*dimVelocity*dimVelocity,
        0.0
    )
);

volScalarField testVegMean_a
    (
        IOobject
        (
            "testVegMean_a",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            ),
        mesh,
        dimensionedScalar
        (
            "testVegMean_a",
            dimensionSet(1, -3, -1, 0, 0, 0, 0),
            0.0
        )
    );

volVectorField UaMean
    (
        IOobject
                (
                    "UaMean",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
            mesh,
            dimensionedVector
                (
                    "UaMean",
                    dimVelocity,
                    vector (0, 0, 0)
                )
        );

///

volScalarField TKEMean_a
(
    IOobject
    (
        "TKEMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
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
        IOobject::NO_READ,
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

volScalarField TKEMeanProdI_a
(
    IOobject
    (
        "TKEMeanProdI_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TKEMeanProdI_a",
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

volScalarField pMean_a
(
    IOobject
    (
        "pMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pMean_a",
        dimPressure,
        0
    )
);

volTensorField UaPrime2MeanF_a
    (
            IOobject
                (
                    "UaPrime2MeanF_a",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
            mesh,
            dimensionedTensor
                (
                    "UaPrime2MeanF_a",
                    dimVelocity*dimVelocity,
                    tensor::zero
                 )
    );

volTensorField alphaSqrUaMean_a
(
    IOobject
    (
        "alphaSqrUaMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "alphaSqrUaMean_a",
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
            IOobject::NO_READ,
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
            IOobject::NO_READ,
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
        IOobject::NO_READ,
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
        IOobject::NO_READ,
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

volScalarField turbDiffusionMeanJ_a
(
    IOobject
    (
        "turbDiffusionMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "turbDiffusionMeanJ_a",
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
        IOobject::NO_READ,
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
        IOobject::NO_READ,
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
        IOobject::NO_READ,
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

volVectorField turbDiffMeanJ_a
(
    IOobject
    (
        "turbDiffMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "turbDiffMeanJ_a",
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
        IOobject::NO_READ,
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

volVectorField viscDiffMeanJ_a
(
    IOobject
    (
        "viscDiffMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "viscDiffMeanJ_a",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volVectorField UaPrimeF_a
    (
            IOobject
                (
                    "UaPrimeF_a",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
            mesh,
            dimensionedVector
                (
                    "UaPrimeF_a",
                    dimVelocity,
                    vector (0, 0, 0)
                )
    );

volVectorField UaMeanF_a
    (
            IOobject
                (
                    "UaMeanF_a",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
            mesh,
            dimensionedVector
                (
                    "UaMeanF_a",
                    dimVelocity,
                    vector (0, 0, 0)
                )
        );

volVectorField alphaUaMean_a
        (
            IOobject
                (
                    "alphaUaMean_a",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
            mesh,
            dimensionedVector
                (
                    "alphaUaMean_a",
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
                    IOobject::NO_READ,
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

volScalarField pPrimeF_a
(
    IOobject
    (
        "pPrimeF_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pPrimeF_a",
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

if (Reset_fluid)
{

Info<<"Manohar in Reset mode = \t"<<endl;

runtime_fluid = 0;

}

else

{

  runtime_fluid = runTime.value();// 0;// this->runTime.timeName();

IOobject betaFluctMean_bHeader
     (
      "betaFluctMean_b",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ
    );

if (betaFluctMean_bHeader.typeHeaderOk<volVectorField>(true))
    {
      volVectorField betaFluctMean_b2
        (
          IOobject
            (
              "betaFluctMean_b",
              runTime.timeName(),
              mesh,
              IOobject::MUST_READ,
              IOobject::NO_WRITE
            ),
          mesh
        );
      betaFluctMean_b = betaFluctMean_b2;
    }

IOobject alpha_aMean_FHeader
    (
      "alpha_aMean_F",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ
    );

if (alpha_aMean_FHeader.typeHeaderOk<volScalarField>(true))
{
  volScalarField alpha_aMean_F2
        (
          IOobject
            (
              "alpha_aMean_F",
              runTime.timeName(),
              mesh,
              IOobject::MUST_READ,
              IOobject::NO_WRITE
            ),
          mesh
        );
        alpha_aMean_F = alpha_aMean_F2;

          IOobject RS2Mean_bHeader
      (
    "RS2Mean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (RS2Mean_bHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField RS2Mean2_b
        (
            IOobject
            (
                "RS2Mean_b",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        RS2Mean_b = RS2Mean2_b;
    }
    // Checks if field RS4Mean is present
    IOobject RS4Mean_bHeader
    (
        "RS4Mean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (RS4Mean_bHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField RS4Mean2_b
        (
            IOobject
            (
                "RS4Mean_b",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        RS4Mean_b = RS4Mean2_b;
    }

    IOobject viscDiffMeanJ_bHeader
    (
        "viscDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (viscDiffMeanJ_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField viscDiffMeanJ2_b
(
    IOobject
    (
        "viscDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
viscDiffMeanJ_b = viscDiffMeanJ2_b;
    }

    IOobject SGSDiffMean_bHeader
    (
        "SGSDiffMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDiffMean_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField SGSDiffMean2_b
(
    IOobject
    (
        "SGSDiffMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDiffMean_b = SGSDiffMean2_b;
    }

   IOobject WbcubeMeanHeader
    (
        "WbcubeMeanMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (WbcubeMeanHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField WbcubeMean2
(
    IOobject
    (
        "WbcubeMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
WbcubeMean = WbcubeMean2;
    }

   IOobject UbsquarWbMeanHeader
    (
        "UbsquarWbMeanMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (UbsquarWbMeanHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField UbsquarWbMean2
(
    IOobject
    (
        "UbsquarWbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
UbsquarWbMean = UbsquarWbMean2;
    }

IOobject turbDiffMeanJ_bHeader
    (
        "turbDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (turbDiffMeanJ_bHeader.typeHeaderOk<volVectorField>(true))
    {

 volVectorField turbDiffMeanJ2_b
(
    IOobject
    (
        "turbDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

turbDiffMeanJ_b = turbDiffMeanJ2_b;
    }

    IOobject dragConstbMeanHeader
    (
        "dragConstbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (dragConstbMeanHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField dragConstbMean2
(
    IOobject
    (
        "dragConstbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

dragConstbMean = dragConstbMean2;

    }

  // Checks if field UaMeanF is present

      // Checks if field UbMeanF_b is present

  IOobject UbMeanF_bHeader
        (
          "UbMeanF_b",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ
        );

  if (UbMeanF_bHeader.typeHeaderOk<volVectorField>(true))
        {
          volVectorField UbMeanF_b2
            (
              IOobject
                (
                  "UbMeanF_b",
                  runTime.timeName(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::NO_WRITE
                ),
              mesh
            );
          UbMeanF_b = UbMeanF_b2;
          betaUbMean_F = UbMeanF_b*(1-alpha_aMean_F);
        }

  IOobject UbMeanHeader
        (
          "UbMean",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ
        );

  if (UbMeanHeader.typeHeaderOk<volVectorField>(true))
        {
          volVectorField UbMean2
            (
              IOobject
                (
                  "UbMean",
                  runTime.timeName(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::NO_WRITE
                ),
              mesh
            );
          UbMean = UbMean2;
        }

      // Checks if field UbPrime2MeanF_b is present

  IOobject UbPrime2MeanF_bHeader
        (
          "UbPrime2MeanF_b",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ
        );

  if (UbPrime2MeanF_bHeader.typeHeaderOk<volTensorField>(true))
        {
          volTensorField UbPrime2MeanF_b2
            (
              IOobject
                (
                  "UbPrime2MeanF_b",
                  runTime.timeName(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::NO_WRITE
                ),
              mesh
            );
          UbPrime2MeanF_b = UbPrime2MeanF_b2;
          betaSqrUbMean_b = (UbPrime2MeanF_b+UbMeanF_b*UbMeanF_b)*
                            (1-alpha_aMean_F);

        }

     // Checks if field UaPrime2MeanF is present

  IOobject viscDissMean_bHeader
    (
      "viscDissMean_b",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ
    );

  if (viscDissMean_bHeader.typeHeaderOk<volScalarField>(true))
    {
      volScalarField viscDissMean2_b
        (
          IOobject
            (
              "viscDissMean_b",
              runTime.timeName(),
              mesh,
              IOobject::MUST_READ,
              IOobject::NO_WRITE
            ),
          mesh
        );
      viscDissMean_b = viscDissMean2_b;
    }

IOobject pMean_bHeader
    (
        "pMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

 if (pMean_bHeader.typeHeaderOk<volScalarField>(true))
    {

    volScalarField pMean2_b
(
    IOobject
    (
        "pMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
  pMean_b = pMean2_b;
    }

     IOobject pPrimeF_bHeader
    (
        "pPrimeF_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pPrimeF_bHeader.typeHeaderOk<volScalarField>(true))
    {

  volScalarField pPrimeF2_b
(
    IOobject
    (
        "pPrimeF_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pPrimeF_b = pPrimeF2_b;
    }

/////////////////

  IOobject viscDiffusionMean_bHeader
    (
        "viscDiffusionMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (viscDiffusionMean_bHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField viscDiffusionMean2_b
(
    IOobject
    (
        "viscDiffusionMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
viscDiffusionMean_b = viscDiffusionMean2_b;
    }

   IOobject SGSDiffusionMean_bHeader
    (
        "SGSDiffusionMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDiffusionMean_bHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField SGSDiffusionMean2_b
(
    IOobject
    (
        "SGSDiffusionMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDiffusionMean_b = SGSDiffusionMean2_b;
    }

     IOobject turbDiffusionMean_bHeader
    (
        "turbDiffusionMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (turbDiffusionMean_bHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField turbDiffusionMean2_b
(
    IOobject
    (
        "turbDiffusionMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
turbDiffusionMean_b = turbDiffusionMean2_b;
    }

IOobject pressDiffusionMean_bHeader
    (
        "pressDiffusionMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pressDiffusionMean_bHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField pressDiffusionMean2_b
(
    IOobject
    (
        "pressDiffusionMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pressDiffusionMean_b = pressDiffusionMean2_b;
    }

     IOobject turbDiffusionMeanDavidHeader
    (
        "turbDiffusionMeanDavid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (turbDiffusionMeanDavidHeader.typeHeaderOk<volVectorField>(true))
    {

volVectorField turbDiffusionMeanDavid2
(
    IOobject
    (
        "turbDiffusionMeanDavid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
turbDiffusionMeanDavid = turbDiffusionMeanDavid2;
    }

/////////////////////////////////

     IOobject SGSDissMean_bHeader
    (
        "SGSDissMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDissMean_bHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField SGSDissMean2_b
(
    IOobject
    (
        "SGSDissMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDissMean_b = SGSDissMean2_b;
    }

///

IOobject beta_bMeanHeader
    (
        "beta_bMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (beta_bMeanHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField beta_bMean2
(
    IOobject
    (
        "beta_bMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
beta_bMean = beta_bMean2;
    }

IOobject TKEMeanProd_bHeader
    (
        "TKEMeanProd_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (TKEMeanProd_bHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField TKEMeanProd_b2
(
    IOobject
    (
        "TKEMeanProd_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
TKEMeanProd_b = TKEMeanProd_b2;
    }

IOobject TKEMean_bHeader
    (
        "TKEMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (TKEMean_bHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField TKEMean_b2
(
    IOobject
    (
        "TKEMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
TKEMean_b = TKEMean_b2;
    }

IOobject alpha_aMean_FHeader
    (
        "alpha_aMean_F",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (alpha_aMean_FHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField alpha_aMean_F2
(
    IOobject
    (
        "alpha_aMean_F",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
alpha_aMean_F = alpha_aMean_F2;
    }

IOobject dragVegMeanHeader
    (
        "dragVegMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (dragVegMeanHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField dragVegMean2
(
    IOobject
    (
        "dragVegMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

dragVegMean = dragVegMean2;

    }

/* Mean Variable for mass budget  */

    IOobject Mean_alphaUaHeader
    (
        "Mean_alphaUa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (Mean_alphaUaHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField Mean_alphaUa2
(
    IOobject
    (
        "Mean_alphaUa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Mean_alphaUa = Mean_alphaUa2;
    }

    IOobject Mean_UaFalphaFHeader
    (
        "Mean_UaFalphaF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (Mean_UaFalphaFHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField Mean_UaFalphaF2
(
    IOobject
    (
        "Mean_UaFalphaF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Mean_UaFalphaF = Mean_UaFalphaF2;
    }

    IOobject MeanUaFHeader
    (
        "MeanUaF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (MeanUaFHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField MeanUaF2
(
    IOobject
    (
        "MeanUaF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
MeanUaF = MeanUaF2;
    }

    IOobject alphaFMeanHeader
        (
            "alphaFMean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
        );

  if (alphaFMeanHeader.typeHeaderOk<volScalarField>(true))
       {
            volScalarField alphaFMean2
              (
                    IOobject
                      (
                          "alphaFMean",
                          runTime.timeName(),
                          mesh,
                          IOobject::MUST_READ,
                          IOobject::NO_WRITE
                      ),
                    mesh
              );
            alphaFMean = alphaFMean2;
       }

 /*   Mean variable defined for Particles             */

          IOobject RS2Mean_aHeader
    (
        "RS2Mean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (RS2Mean_aHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField RS2Mean2_a
        (
            IOobject
            (
                "RS2Mean_a",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        RS2Mean_a = RS2Mean2_a;
    }
    // Checks if field RS4Mean is present
    IOobject RS4Mean_aHeader
    (
        "RS4Mean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (RS4Mean_aHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField RS4Mean2_a
        (
            IOobject
            (
                "RS4Mean_a",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        RS4Mean_a = RS4Mean2_a;
    }

    IOobject viscDiffMeanJ_aHeader
    (
        "viscDiffMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (viscDiffMeanJ_aHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField viscDiffMeanJ2_a
(
    IOobject
    (
        "viscDiffMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
viscDiffMeanJ_a = viscDiffMeanJ2_a;
    }

     IOobject turbDiffMean_aHeader
    (
        "turbDiffMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (turbDiffMean_aHeader.typeHeaderOk<volVectorField>(true))
    {

 volVectorField turbDiffMean2_a
(
    IOobject
    (
        "turbDiffMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
turbDiffMean_a = turbDiffMean2_a;
    }

IOobject turbDiffMeanJ_aHeader
    (
        "turbDiffMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (turbDiffMeanJ_aHeader.typeHeaderOk<volVectorField>(true))
    {

 volVectorField turbDiffMeanJ2_a
(
    IOobject
    (
        "turbDiffMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

turbDiffMeanJ_a = turbDiffMeanJ2_a;
    }

    // Checks if field UaMeanF_a is present

            // Checks if field UbMeanF is present

    IOobject UaMeanF_aHeader
                    (
                        "UaMeanF_a",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    );

    if (UaMeanF_aHeader.typeHeaderOk<volVectorField>(true))
                    {
                        volVectorField UaMeanF_a2
                            (
                                    IOobject
                                        (
                                            "UaMeanF_a",
                                            runTime.timeName(),
                                            mesh,
                                            IOobject::MUST_READ,
                                            IOobject::NO_WRITE
                                        ),
                                    mesh
                            );
                        UaMeanF_a = UaMeanF_a2;
                        alphaUaMean_a = UaMeanF_a*(alpha_aMean_a);
                    }

            // Checks if field UbPrime2MeanF is present

    IOobject UaPrime2MeanF_aHeader
                    (
                        "UaPrime2MeanF_a",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    );

    if (UaPrime2MeanF_aHeader.typeHeaderOk<volTensorField>(true))
                    {
                        volTensorField UaPrime2MeanF_a2
                            (
                                    IOobject
                                        (
                                            "UaPrime2MeanF_a",
                                            runTime.timeName(),
                                            mesh,
                                            IOobject::MUST_READ,
                                            IOobject::NO_WRITE
                                        ),
                                    mesh
                            );
                       UaPrime2MeanF_a = UaPrime2MeanF_a2;
                       alphaSqrUaMean_a = (UaPrime2MeanF_a+UaMeanF_a*UaMeanF_a)
                                          *(alpha_aMean_a);

                    }

       // Checks if field UaPrime2MeanF_a is present

    IOobject viscDissMean_aHeader
            (
                "viscDissMean_a",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            );

    if (viscDissMean_aHeader.typeHeaderOk<volScalarField>(true))
            {
                volScalarField viscDissMean2_a
                    (
                            IOobject
                                (
                                    "viscDissMean_a",
                                    runTime.timeName(),
                                    mesh,
                                    IOobject::MUST_READ,
                                    IOobject::NO_WRITE
                                ),
                            mesh
                    );
                viscDissMean_a = viscDissMean2_a;
            }

       // Checks if field UaPrime2MeanF_a is present

IOobject pMean_aHeader
    (
        "pMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

 if (pMean_aHeader.typeHeaderOk<volScalarField>(true))
    {

    volScalarField pMean2_a
(
    IOobject
    (
        "pMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
  pMean_a = pMean2_a;
    }

     IOobject pPrimeF_aHeader
    (
        "pPrimeF_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pPrimeF_aHeader.typeHeaderOk<volScalarField>(true))
    {

  volScalarField pPrimeF2_a
(
    IOobject
    (
        "pPrimeF_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pPrimeF_a = pPrimeF2_a;
    }

     IOobject SGSDissMean_aHeader
    (
        "SGSDissMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDissMean_aHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField SGSDissMean2_a
(
    IOobject
    (
        "SGSDissMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDissMean_a = SGSDissMean2_a;
    }

/////////////////

  IOobject viscDiffusionMean_aHeader
    (
        "viscDiffusionMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (viscDiffusionMean_aHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField viscDiffusionMean2_a
(
    IOobject
    (
        "viscDiffusionMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
viscDiffusionMean_a = viscDiffusionMean2_a;
    }

   IOobject SGSDiffusionMean_aHeader
    (
        "SGSDiffusionMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDiffusionMean_aHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField SGSDiffusionMean2_a
(
    IOobject
    (
        "SGSDiffusionMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDiffusionMean_a = SGSDiffusionMean2_a;
    }

     IOobject turbDiffusionMean_aHeader
    (
        "turbDiffusionMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (turbDiffusionMean_aHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField turbDiffusionMean2_a
(
    IOobject
    (
        "turbDiffusionMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
turbDiffusionMean_a = turbDiffusionMean2_a;
    }

IOobject pressDiffusionMean_aHeader
    (
        "pressDiffusionMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pressDiffusionMean_aHeader.typeHeaderOk<volScalarField>(true))
    {

volScalarField pressDiffusionMean2_a
(
    IOobject
    (
        "pressDiffusionMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pressDiffusionMean_a = pressDiffusionMean2_a;
    }

  IOobject pMech_bHeader
    (
        "pMech_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pMech_bHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField pMech2
(
    IOobject
    (
        "pMech_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pMech_b = pMech2;
    }

IOobject pAHeader
    (
        "pA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pAHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField pA2
(
    IOobject
    (
        "pA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pA = pA2;
    }

IOobject pFFHeader
    (
        "pFF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pFFHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField pFF2
(
    IOobject
    (
        "pFF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pFF = pFF2;
    }

IOobject rho_mixdHeader
    (
        "rho_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (rho_mixdHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField rho_mixd2
(
    IOobject
    (
        "rho_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
rho_mixd = rho_mixd2;
    }

IOobject nonl_mixdHeader
    (
        "nonl_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (nonl_mixdHeader.typeHeaderOk<volTensorField>(true))
    {
volTensorField nonl_mixd2
(
    IOobject
    (
        "nonl_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
nonl_mixd = nonl_mixd2;
    }

IOobject buoy_mixdHeader
    (
        "buoy_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (buoy_mixdHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField buoy_mixd2
(
    IOobject
    (
        "buoy_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
buoy_mixd = buoy_mixd2;
    }

IOobject SolidPressureContr_bHeader
    (
        "SolidPressureContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SolidPressureContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField SolidPressureContr2
(
    IOobject
    (
        "SolidPressureContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SolidPressureContr_b = SolidPressureContr2;
    }

IOobject divTauaContr_bHeader
    (
        "divTauaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (divTauaContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField divTauaContr2
(
    IOobject
    (
        "divTauaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
divTauaContr_b = divTauaContr2;
    }

IOobject ViscStressContr_bHeader
    (
        "ViscStressContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (ViscStressContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField ViscStressContr2
(
    IOobject
    (
        "ViscStressContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
ViscStressContr_b = ViscStressContr2;
    }

IOobject FrictContr_bHeader
    (
        "FrictContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (FrictContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField FrictContr2
(
    IOobject
    (
        "FrictContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
FrictContr_b = FrictContr2;
    }

IOobject phiRaContr_bHeader
    (
        "phiRaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (phiRaContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField phiRaContr2
(
    IOobject
    (
        "phiRaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
phiRaContr_b = phiRaContr2;
    }

IOobject phiRa2Contr_bHeader
    (
        "phiRa2Contr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (phiRa2Contr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField phiRa2Contr2
(
    IOobject
    (
        "phiRa2Contr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
phiRa2Contr_b = phiRa2Contr2;
    }

IOobject divRcaContr_bHeader
    (
        "divRcaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (divRcaContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField divRcaContr2
(
    IOobject
    (
        "divRcaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
divRcaContr_b = divRcaContr2;
    }

IOobject RcaContr_bHeader
    (
        "RcaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (RcaContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField RcaContr2
(
    IOobject
    (
        "RcaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
RcaContr_b = RcaContr2;
    }

IOobject BulkContr_bHeader
    (
        "BulkContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (BulkContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField BulkContr2
(
    IOobject
    (
        "BulkContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
BulkContr_b = BulkContr2;
    }

IOobject gravityContr_bHeader
    (
        "gravityContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (gravityContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField gravityContr2
(
    IOobject
    (
        "gravityContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
gravityContr_b = gravityContr2;
    }

IOobject DragContr_bHeader
    (
        "DragContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (DragContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField DragContr2
(
    IOobject
    (
        "DragContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
DragContr_b = DragContr2;
    }

IOobject SuspContr_bHeader
    (
        "SuspContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SuspContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField SuspContr2
(
    IOobject
    (
        "SuspContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SuspContr_b = SuspContr2;
    }

IOobject ExtContr_bHeader
    (
        "ExtContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (ExtContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField ExtContr2
(
    IOobject
    (
        "ExtContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
ExtContr_b = ExtContr2;
    }

IOobject ErosionHeader
    (
        "Erosion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (ErosionHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField Erosion2
(
    IOobject
    (
        "Erosion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Erosion = Erosion2;
    }

IOobject DepositionHeader
    (
        "Deposition",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (DepositionHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField Deposition2
(
    IOobject
    (
        "Deposition",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Deposition = Deposition2;
    }

/////////////////

IOobject TauaHeader
    (
        "Taua",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (TauaHeader.typeHeaderOk<volTensorField>(true))
    {
volTensorField Taua2
(
    IOobject
    (
        "Taua",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Taua = Taua2;
    }

IOobject TaubHeader
    (
        "Taub",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (TaubHeader.typeHeaderOk<volTensorField>(true))
    {
volTensorField Taub2
(
    IOobject
    (
        "Taub",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Taub = Taub2;
    }

}

}
