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
    format      binary;
    class       volScalarField;
    location    "0";
    object      k;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0; 

boundaryField
{
    cylinder
    {
        type            kqRWallFunction;
        value           uniform 1e-10; 
    }
    inlet
    {
        type            groovyBC;
        refValue        uniform 0;
        value           uniform 1e-10;
        refGradient     uniform 0;
        valueFraction   uniform 1;
        valueExpression "inletprofilek(zp)";
        gradientExpression "0";
        fractionExpression "1";
        evaluateDuringConstruction 1;
        variables       "zp=pos().z;";
        timelines       (
);
        lookuptables    (
{
        name            inletprofilek;
        file        "$FOAM_CASE/1d_profil/k.b.xy";
        outOfBounds     clamp;
}
);
    }
    outlet
    {
        type            zeroGradient;
    }
    outletb
    {
        type            zeroGradient;
    }
    pit
    {
        type            fixedValue;
        value           uniform 1e-12;
    }
    lateral
    {
        type            cyclic;
    }
    symplane
    {
        type            cyclic;
    }
    bottom
    {
        type            zeroGradient;
    }
    surface
    {
        type            zeroGradient;
    }
}


// ************************************************************************* //
