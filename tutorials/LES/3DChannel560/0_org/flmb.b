/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5                                     |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "1";
    object      flm;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 4 -4 0 0 0 0];

internalField   uniform 0.1;

boundaryField
{
    bottomWall
    {
        type            fixedValue;
        value           uniform 0;
    }
    topWall
    {
        type            fixedValue;
        value           uniform 0;
    }
    sides1_half0
    {
        type            cyclic;
    }
    sides2_half0
    {
        type            cyclic;
    }
    inout1_half0
    {
        type            cyclic;
    }
    inout2_half0
    {
        type            cyclic;
    }
    sides2_half1
    {
        type            cyclic;
    }
    sides1_half1
    {
        type            cyclic;
    }
    inout1_half1
    {
        type            cyclic;
    }
    inout2_half1
    {
        type            cyclic;
    }
}


// ************************************************************************* //
