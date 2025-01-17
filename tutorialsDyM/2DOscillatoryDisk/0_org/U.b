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

    walls
    {
        //type            slip;
     type            uniformFixedValue;
        value           uniform (0 0 0);    
    }

    hole
    {
        type            movingWallVelocity;
        value           uniform (0 0 0);
    }

//    left1
//    {
//        type            pressureInletOutletVelocity;
//        value           uniform (0 0 0);
//    }
//    left1
//    {
//        type            fixedValue;
//        value           $internalField;
//    }
    top
    {
     type            zeroGradient;

    }
    bottom
    {
     type            fixedValue;
        value           uniform (0 0 0);    }
    inlet
    {
     type            fixedValue;
        value           uniform (0 0 0);    }    
    outlet
    {
     type            fixedValue;
        value           uniform (0 0 0);    
    }

    overset
    {
        type            overset;
    }
}

// ************************************************************************* //
