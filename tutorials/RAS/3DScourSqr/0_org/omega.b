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
    object      omega;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

internalField   uniform 1882; 

boundaryField
{
    cylinder
    {
        type            omegaWallFunction;
	blending	tanh;
	value           uniform 0.033;
    }
    inlet
    {
        type            groovyBC;
        refvalue           uniform 0.033;
        value           uniform 0.033;
        refGradient     uniform 0;
        valueFraction   uniform 1;
        valueExpression "inletprofileomega(pos().z)";
        gradientExpression "0";
        fractionExpression "1";
        evaluateDuringConstruction 1;
        variables       "";
        timelines       (
);
        lookuptables    (
{
        name            inletprofileomega;
        file        "$FOAM_CASE/1d_profil/omega.b.xy";
        outOfBounds     clamp;
}
);
    }
    outlet
    {
        type            zeroGradient;
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
        type            fixedValue;
        value           uniform 100;
    }
    surface
    {
        type            zeroGradient;
    }
}


// ************************************************************************* //
