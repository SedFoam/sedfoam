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
/*
dimensionedScalar zeroo
(
	"zeroo",
	dimensionSet(0, 2, -2, 0, 0, 0, 0),
	scalar(0.0)
);

*/

Switch favreAverage_particle
(
    runTime.controlDict().lookupOrDefault("favreAverage_particle", false)
);


Switch Reset_particle
(
	runTime.controlDict().lookupOrDefault("Reset_particle", false)
);


if (Reset_particle)
{
	Info<<"Reset particle options = \t"<<Reset_particle<<endl;
}


Info<<"Reset particle Mode = \t"<<Reset_particle<<endl;

scalar timeaverage_particle = 0;
scalar runtime_particle     = 0;
double time_diff_particle = 0;
double time_integrate_particle = 0;

//scalar deltaT = 0;

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

volScalarField SGSDissMeanJ_a
(
    IOobject
    (
        "SGSDissMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDissMeanJ_a",
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


volScalarField pressDiffusionMeankin_a
(
    IOobject
    (
        "pressDiffusionMeankin_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pressDiffusionMeankin_a",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);


volScalarField pressDiffusionMeancon_a
(
    IOobject
    (
        "pressDiffusionMeancon_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pressDiffusionMeancon_a",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField pressDiffusionMeanJ_a
(
    IOobject
    (
        "pressDiffusionMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pressDiffusionMeanJ_a",
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

volScalarField SGSDiffusionMeanJ_a
(
    IOobject
    (
        "SGSDiffusionMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDiffusionMeanJ_a",
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

volVectorField pressDiffMean_a
(
    IOobject
    (
        "pressDiffMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "pressDiffMean_a",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);


volVectorField pressDiffMeankin_a
(
    IOobject
    (
        "pressDiffMeankin_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "pressDiffMeankin_a",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);


volVectorField pressDiffMeancon_a
(
    IOobject
    (
        "pressDiffMeancon_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "pressDiffMeancon_a",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volVectorField pressDiffMeanJ_a
(
    IOobject
    (
        "pressDiffMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "pressDiffMeanJ_a",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
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

volScalarField viscDiffusionMeanJ_a
(
    IOobject
    (
        "viscDiffusionMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "viscDiffusionMeanJ_a",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volVectorField SGSDiffMean_a
(
    IOobject
    (
        "SGSDiffMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "SGSDiffMean_a",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volVectorField SGSDiffMeanJ_a
(
    IOobject
    (
        "SGSDiffMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "SGSDiffMeanJ_a",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volVectorField viscDiffMean_a
(
    IOobject
    (
        "viscDiffMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "viscDiffMean_a",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
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

//~ volScalarField viscDissMeankin_a
	//~ (
   			//~ IOobject
   				//~ (
       				//~ "viscDissMeankin_a",
       				//~ runTime.timeName(),
       				//~ mesh,
       				//~ IOobject::NO_READ,
       				//~ IOobject::NO_WRITE
   				//~ ),
   			//~ mesh,
   			//~ dimensionedScalar
   				//~ (
       				//~ "viscDissMeankin_a",
       				//~ dimVelocity*dimVelocity*dimVelocity/dimLength,
       				//~ 0.0
   				//~ )
	//~ );

volScalarField viscDissMeancon_a
	(
   			IOobject
   				(
       				"viscDissMeancon_a",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
   				),
   			mesh,
   			dimensionedScalar
   				(
       				"viscDissMeancon_a",
       				dimVelocity*dimVelocity*dimVelocity/dimLength,
       				0.0
   				)
	);


volScalarField viscDissMeanJ_a
	(
   			IOobject
   				(
       				"viscDissMeanJ_a",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
   				),
   			mesh,
   			dimensionedScalar
   				(
       				"viscDissMeanJ_a",
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



if (Reset_particle)
{

Info<<"Manohar in Reset mode = \t"<<endl;

runtime_particle = 0;

}

else 

{

	runtime_particle = runTime.value();// 0;// this->runTime.timeName();

	Info<<"runtime_particle in create file =\t"<<runtime_particle<<endl;
	Info<<"Manohar is not in Reset Mode = \t"<<runTime.value()<<endl;


IOobject alpha_aMeanHeader
    		(
        		"alpha_aMean",
        		runTime.timeName(),
        		mesh,
        		IOobject::NO_READ
    		);


if (alpha_aMeanHeader.typeHeaderOk<volScalarField>(true))
{
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

    IOobject viscDiffMean_aHeader
    (
        "viscDiffMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (viscDiffMean_aHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField viscDiffMean2_a
(
    IOobject
    (
        "viscDiffMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
viscDiffMean_a = viscDiffMean2_a;
}


    IOobject viscDiffMeanJ_aHeader
    (
        "viscDiffMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (viscDiffMeanJ_aHeader.typeHeaderOk<volTensorField>(true))
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


    IOobject SGSDiffMean_aHeader
    (
        "SGSDiffMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDiffMean_aHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField SGSDiffMean2_a
(
    IOobject
    (
        "SGSDiffMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDiffMean_a = SGSDiffMean2_a;
}


    IOobject SGSDiffMeanJ_aHeader
    (
        "SGSDiffMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDiffMeanJ_aHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField SGSDiffMeanJ2_a
(
    IOobject
    (
        "SGSDiffMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDiffMeanJ_a = SGSDiffMeanJ2_a;
}

     IOobject turbDiffMean_aHeader
    (
        "turbDiffMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (turbDiffMean_aHeader.typeHeaderOk<volTensorField>(true))
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
    if (turbDiffMeanJ_aHeader.typeHeaderOk<volTensorField>(true))
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

     IOobject pressDiffMean_aHeader
    (
        "pressDiffMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pressDiffMean_aHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField pressDiffMean2_a
(
    IOobject
    (
        "pressDiffMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pressDiffMean_a = pressDiffMean2_a;
}



     IOobject pressDiffMeankin_aHeader
    (
        "pressDiffMeankin_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pressDiffMeankin_aHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField pressDiffMeankin2_a
(
    IOobject
    (
        "pressDiffMeankin_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pressDiffMeankin_a = pressDiffMeankin2_a;
}


     IOobject pressDiffMeancon_aHeader
    (
        "pressDiffMeancon_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pressDiffMeancon_aHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField pressDiffMeancon2_a
(
    IOobject
    (
        "pressDiffMeancon_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pressDiffMeancon_a = pressDiffMeancon2_a;
}

     IOobject pressDiffMeanJ_aHeader
    (
        "pressDiffMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pressDiffMeanJ_aHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField pressDiffMeanJ2_a
(
    IOobject
    (
        "pressDiffMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pressDiffMeanJ_a = pressDiffMeanJ2_a;
}


	// Checks if field UaMeanF is present
	
			// Checks if field UbMeanF is present
	
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
        				alphaUaMean = UaMeanF*(alpha_aMean);
    				}

		
	    
	
			// Checks if field UbPrime2MeanF is present
	
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
				        alphaSqrUaMean = (UaPrime2MeanF+UaMeanF*UaMeanF)*(alpha_aMean);


    				}

	   // Checks if field UaPrime2MeanF is present
    				
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
    		
	   // Checks if field UaPrime2MeanF is present
    				
  	//~ IOobject viscDissMeankin_aHeader
    		//~ (
        		//~ "viscDissMeankin_a",
        		//~ runTime.timeName(),
        		//~ mesh,
        		//~ IOobject::NO_READ
    		//~ );
    	
	//~ if (viscDissMeankin_aHeader.typeHeaderOk<volScalarField>(true))
    		//~ {
        		//~ volScalarField viscDissMeankin2_a
        			//~ (
            				//~ IOobject
            					//~ (
                					//~ "viscDissMeankin_a",
                					//~ runTime.timeName(),
                					//~ mesh,
                					//~ IOobject::MUST_READ,
                					//~ IOobject::NO_WRITE
            					//~ ),
            				//~ mesh
        			//~ );
        		//~ viscDissMeankin_a = viscDissMeankin2_a;
    		//~ }

	   // Checks if field UaPrime2MeanF is present
    				
  	IOobject viscDissMeancon_aHeader
    		(
        		"viscDissMeancon_a",
        		runTime.timeName(),
        		mesh,
        		IOobject::NO_READ
    		);
    	
	if (viscDissMeancon_aHeader.typeHeaderOk<volScalarField>(true))
    		{
        		volScalarField viscDissMeancon2_a
        			(
            				IOobject
            					(
                					"viscDissMeancon_a",
                					runTime.timeName(),
                					mesh,
                					IOobject::MUST_READ,
                					IOobject::NO_WRITE
            					),
            				mesh
        			);
        		viscDissMeancon_a = viscDissMeancon2_a;
    		}
    		
  	IOobject viscDissMeanJ_aHeader
    		(
        		"viscDissMeanJ_a",
        		runTime.timeName(),
        		mesh,
        		IOobject::NO_READ
    		);
    	
	if (viscDissMeanJ_aHeader.typeHeaderOk<volScalarField>(true))
    		{
        		volScalarField viscDissMeanJ2_a
        			(
            				IOobject
            					(
                					"viscDissMeanJ_a",
                					runTime.timeName(),
                					mesh,
                					IOobject::MUST_READ,
                					IOobject::NO_WRITE
            					),
            				mesh
        			);
        		viscDissMeanJ_a = viscDissMeanJ2_a;
    		}



IOobject pMean_aHeader
    (
        "pMean_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

 if (pMean_aHeader.typeHeaderOk<volTensorField>(true))
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
    if (pPrimeF_aHeader.typeHeaderOk<volTensorField>(true))
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
    if (SGSDissMean_aHeader.typeHeaderOk<volTensorField>(true))
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


     IOobject SGSDissMeanJ_aHeader
    (
        "SGSDissMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDissMeanJ_aHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField SGSDissMeanJ2_a
(
    IOobject
    (
        "SGSDissMeanJ_a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDissMeanJ_a = SGSDissMeanJ2_a;
}

}

}


