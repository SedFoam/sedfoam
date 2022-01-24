/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5.0                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      binary;
    class       volScalarField;
    location    "0";
    object      nut.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0; 

boundaryField
{
    bottomWall
    {
        type            zeroGradient;
    }
    topWall
    {
        type            zeroGradient;
    }
    sides1
    {
        type            cyclic;
    }
    sides2
    {
        type            cyclic;
    }
    inlet
    {
        type            cyclic;
    }
    outlet
    {
        type            cyclic;
    }
}


// ************************************************************************* //