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
    object      nut;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 1e-10;

boundaryField
{
    cylinder
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
    inlet
    {
        type            cyclic;
    }
    outlet
    {
        type            cyclic;
    }
    surface
    {
        type            zeroGradient;
    }
    bottom
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
    lateralfront
    {
        type            empty;
    }
    lateralback
    {
        type            empty;
    }
}


// ************************************************************************* //
