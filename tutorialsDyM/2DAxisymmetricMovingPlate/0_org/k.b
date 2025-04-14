/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1906                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      delta;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 1e-6;

boundaryField
{
    front
    {
         type            wedge;

    }
    back
    {
         type            wedge;

    }

    hole
    {
        type		zeroGradient;
    }
    top
    {
//        type            fixedValue;
  //              value uniform 0;
        type            zeroGradient;

    }
    bottom
    {
        type            zeroGradient;
    }
    
    
    "overset.*"
    {
        type            overset;
    }
    OuterWall
    {
        type            zeroGradient;
    }
    walls
    {
        type            zeroGradient;
    }
    axis
    {
        type            empty;
    }
}


// ************************************************************************* //
