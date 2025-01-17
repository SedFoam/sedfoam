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
    class       volScalarField;
    location    "0";
    object      k;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField     uniform 1e-6;

boundaryField
{
    cylinder
    {
        type            kqRWallFunction;
        value           uniform 1e-10;
    }
    inlet
    {
        type            cyclic;
    }
    outlet
    {
        type            cyclic;
    }
    lateralfront
    {
        type            empty;
    }
    lateralback
    {
        type            empty;
    }
    bottom
    {
        type            kqRWallFunction;
        value           uniform 1e-10;
    }
    surface
    {
        type            zeroGradient;
    }
}


// ************************************************************************* //
