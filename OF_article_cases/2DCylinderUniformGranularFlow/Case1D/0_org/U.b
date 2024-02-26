/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.4.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "0";
    object      Ua;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (1.e-2 0 0);

boundaryField
{
    cylinder
    {

           type            slip;

//        type            fixedValue;
  //      value           uniform (0 0 0);
    }
    inlet
    {
        type            cyclic;

		
    }
    outlet
    {
        type            cyclic;

    }
    top
    {
      //  type            fixedValue;
       // value           uniform (0.01 0 0);
//		type            slip;
type zeroGradient;  
  }
    bottom
    {
        type            fixedValue;
        value           uniform (0.01 0 0);
//		type            slip;

    }
    frontAndBack
    {
        type            empty;
    }
}


// ************************************************************************* //
