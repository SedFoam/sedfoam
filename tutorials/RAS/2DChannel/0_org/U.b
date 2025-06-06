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

internalField   uniform (0.01 0 0);

boundaryField
{
    inlet
    {
        type            groovyBC;
        refValue        uniform  (0 0 0);
        value           uniform  (0 0 0);
        variables       "t=time();";
	valueExpression "max(min(t/10.,1),0.01)*vector(1,0,0)";
        gradientExpression "vector(0,0,0)";
        fractionExpression "1";
    }
    outlet
    {
      type            pressureInletOutletVelocity;
        value           uniform (0 0 0);
    }
    top
    {
        type            slip;
    }
    walls
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
    apron
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
    front
    {
        type            empty;
    }
    back
    {
        type            empty;
    }
}


// ************************************************************************* //
