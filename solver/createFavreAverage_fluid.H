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
double runtime_sub     = 0;
scalar deltaT = 0;
scalar T0 = 0;
double runtime_Manohar = 0;
double time_integrate_fluid = 0;

double total_time = 0;
double time_diff_fluid = 0;

if (AverageMode)
{
	T0 =runTime.controlDict().get<scalar>("StartAverageTime");
}

Info<<"To before else looop in Average fluid non reset mode = \t"<<T0<<endl;

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


volScalarField turbDiffMeanDavid
(
    IOobject
    (
        "turbDiffMeanDavid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "turbDiffMeanDavid",
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

volScalarField SGSDissMeanJ_b
(
    IOobject
    (
        "SGSDissMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDissMeanJ_b",
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

volScalarField turbDiffusionMeanJ_b
(
    IOobject
    (
        "turbDiffusionMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "turbDiffusionMeanJ_b",
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

volScalarField pressDiffusionMeanJ_b
(
    IOobject
    (
        "pressDiffusionMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pressDiffusionMeanJ_b",
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

volScalarField SGSDiffusionMeanJ_b
(
    IOobject
    (
        "SGSDiffusionMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDiffusionMeanJ_b",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volVectorField turbDiffMean_b
(
    IOobject
    (
        "turbDiffMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "turbDiffMean_b",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
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

volVectorField pressDiffMean_b
(
    IOobject
    (
        "pressDiffMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "pressDiffMean_b",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);


volVectorField pressDiffMeanJ_b
(
    IOobject
    (
        "pressDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "pressDiffMeanJ_b",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
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

volScalarField viscDiffusionMeanJ_b
(
    IOobject
    (
        "viscDiffusionMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "viscDiffusionMeanJ_b",
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

volVectorField SGSDiffMeanJ_b
(
    IOobject
    (
        "SGSDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "SGSDiffMeanJ_b",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volVectorField viscDiffMean_b
(
    IOobject
    (
        "viscDiffMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "viscDiffMean_b",
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

volScalarField viscDissMeanJ_b
	(
   			IOobject
   				(
       				"viscDissMeanJ_b",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
   				),
   			mesh,
   			dimensionedScalar
   				(
       				"viscDissMeanJ_b",
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

    IOobject viscDiffMean_bHeader
    (
        "viscDiffMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (viscDiffMean_bHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField viscDiffMean2_b
(
    IOobject
    (
        "viscDiffMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
viscDiffMean_b = viscDiffMean2_b;
}


    IOobject viscDiffMeanJ_bHeader
    (
        "viscDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (viscDiffMeanJ_bHeader.typeHeaderOk<volTensorField>(true))
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
    if (SGSDiffMean_bHeader.typeHeaderOk<volTensorField>(true))
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


    IOobject SGSDiffMeanJ_bHeader
    (
        "SGSDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDiffMeanJ_bHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField SGSDiffMeanJ2_b
(
    IOobject
    (
        "SGSDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDiffMeanJ_b = SGSDiffMeanJ2_b;
}


   IOobject WbcubeMeanHeader
    (
        "WbcubeMeanMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (WbcubeMeanHeader.typeHeaderOk<volTensorField>(true))
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
    if (UbsquarWbMeanHeader.typeHeaderOk<volTensorField>(true))
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



     IOobject turbDiffMean_bHeader
    (
        "turbDiffMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (turbDiffMean_bHeader.typeHeaderOk<volTensorField>(true))
    {
 
 volVectorField turbDiffMean2_b
(
    IOobject
    (
        "turbDiffMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
turbDiffMean_b = turbDiffMean2_b;
}


IOobject turbDiffMeanJ_bHeader
    (
        "turbDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (turbDiffMeanJ_bHeader.typeHeaderOk<volTensorField>(true))
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



 IOobject turbDiffMeanDavidHeader
    (
        "turbDiffMeanDavid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (turbDiffMeanDavidHeader.typeHeaderOk<volTensorField>(true))
    {
 
 volScalarField turbDiffMean2David
(
    IOobject
    (
        "turbDiffMeanDavid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
turbDiffMeanDavid = turbDiffMean2David;
}


    IOobject dragConstbMeanHeader
    (
        "dragConstbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (dragConstbMeanHeader.typeHeaderOk<volTensorField>(true))
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

IOobject pressDiffMean_bHeader
    (
        "pressDiffMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pressDiffMean_bHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField pressDiffMean2_b
(
    IOobject
    (
        "pressDiffMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pressDiffMean_b = pressDiffMean2_b;
}


     IOobject pressDiffMeanJ_bHeader
    (
        "pressDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pressDiffMeanJ_bHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField pressDiffMeanJ2_b
(
    IOobject
    (
        "pressDiffMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pressDiffMeanJ_b = pressDiffMeanJ2_b;
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
				        betaSqrUbMean_b = (UbPrime2MeanF_b+UbMeanF_b*UbMeanF_b)*(1-alpha_aMean_F);


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
    		
    		
  	IOobject viscDissMeanJ_bHeader
    		(
        		"viscDissMeanJ_b",
        		runTime.timeName(),
        		mesh,
        		IOobject::NO_READ
    		);
    	
	if (viscDissMeanJ_bHeader.typeHeaderOk<volScalarField>(true))
    		{
        		volScalarField viscDissMeanJ2_b
        			(
            				IOobject
            					(
                					"viscDissMeanJ_b",
                					runTime.timeName(),
                					mesh,
                					IOobject::MUST_READ,
                					IOobject::NO_WRITE
            					),
            				mesh
        			);
        		viscDissMeanJ_b = viscDissMeanJ2_b;
    		}



IOobject pMean_bHeader
    (
        "pMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

 if (pMean_bHeader.typeHeaderOk<volTensorField>(true))
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
    if (pPrimeF_bHeader.typeHeaderOk<volTensorField>(true))
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


     IOobject SGSDissMean_bHeader
    (
        "SGSDissMean_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDissMean_bHeader.typeHeaderOk<volTensorField>(true))
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


     IOobject SGSDissMeanJ_bHeader
    (
        "SGSDissMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDissMeanJ_bHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField SGSDissMeanJ2_b
(
    IOobject
    (
        "SGSDissMeanJ_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDissMeanJ_b = SGSDissMeanJ2_b;
}

///

IOobject beta_bMeanHeader
    (
        "beta_bMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (beta_bMeanHeader.typeHeaderOk<volTensorField>(true))
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

IOobject alpha_aMean_FHeader
    (
        "alpha_aMean_F",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (alpha_aMean_FHeader.typeHeaderOk<volTensorField>(true))
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
			
}

}


