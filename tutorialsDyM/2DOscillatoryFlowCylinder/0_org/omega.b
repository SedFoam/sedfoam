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
    object      omega;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

internalField   uniform 10.; 

boundaryField
{
    cylinder
    {
        type            omegaWallFunction;
        value           uniform 10.;
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
    surface
    {
        type            zeroGradient;
    }
    bottom
    {
        type            omegaWallFunction;
        value           uniform 10.;
    }
    defaultFaces
    {
        type            empty;
    }
}


// ************************************************************************* //
