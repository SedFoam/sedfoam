/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2006                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{
    #includeEtc "caseDicts/setConstraintTypes"
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
        type            movingWallVelocity;
        value           uniform (0 0 0);    }
    top
    {
		
        type            zeroGradient;

    }
    bottom
    {
     type            fixedValue;
        value           uniform (0 0 0);    }
    
    overset
    {
        type            overset;
    }
    OuterWall
    {
     type            fixedValue;
        value           uniform (0 0 0);    }
    walls
    {
     type            uniformFixedValue;
        value           uniform (0 0 0);    }
    axis
    {
        type            empty;
    }
    
   
}

// ************************************************************************* //
