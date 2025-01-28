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

Switch favreAverage_mass
(
    runTime.controlDict().lookupOrDefault("favreAverage_mass", false)
);


Switch Reset_mass
(
	runTime.controlDict().lookupOrDefault("Reset_mass", false)
);


if (Reset_mass)
{
	Info<<"Reset Mass options = \t"<<Reset_mass<<endl;
}


Info<<"Reset Mass Mode = \t"<<Reset_mass<<endl;

//scalar timeaverage = 0;
//scalar runtime     = 0;
//scalar deltaT = 0;

scalar timeaverage_mass = 0;
scalar runtime_mass     = 0;
double time_diff_mass = 0;
double time_integrate_mass = 0;

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

if (Reset_particle)
{

Info<<"Manohar in Reset mode = \t"<<endl;

runtime_mass = 0;

}

else 

{

	runtime_mass = runTime.value();// 0;// this->runTime.timeName();

	Info<<"runtime_mass in create file =\t"<<runtime_mass<<endl;
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

		
///
			
}

}


