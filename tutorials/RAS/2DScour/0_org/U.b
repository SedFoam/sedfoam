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
    object      Ub;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{
    inandouthalf11
    {
        type            groovyBC;
        refValue        uniform ( 0 0 0);
        refGradient     uniform (0 0 0);
        valueFraction   uniform 1;
        value           uniform (0 0 0);
         variables "yp=pos().y;";
        valueExpression "(0.0369/0.41)*log(30.0*mag(yp)/(2.5*0.00025))*vector(1,0,0)";
        gradientExpression "vector(0,0,0)";
        fractionExpression "1";
        evaluateDuringConstruction 1;
    }
    inandouthalf12
    {
        type            inletOutlet;
        inletValue      uniform (0 0 0);
        value           uniform (0 0 0);
    }
    inandouthalf21
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
    inandouthalf22
    {
      type            fixedValue;
      value           uniform (0 0 0);
    }
    top
    {
        type            zeroGradient;
    }
    walls
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
    frontAndBackPlanes
    {
        type            empty;
    }
}


// ************************************************************************* //
