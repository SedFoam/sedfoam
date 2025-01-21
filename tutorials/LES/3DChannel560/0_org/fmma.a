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
    object      fmm;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 4 -4 0 0 0 0];

internalField   uniform 0.1;

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
    "(left|right)"
    {
        type            cyclic;
    }
    "(inlet|outlet)"
    {
        type            cyclic;
    }
}


// ************************************************************************* //
