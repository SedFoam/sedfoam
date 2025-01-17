/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2306                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "0";
    object      Ua;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (1.e-2 0 0);

boundaryField
{
    cylinder
    {
   type            movingWallVelocity;
   value           uniform (0 0 0);
    }
    inlet
    {
        type            cyclic;
    }
    outlet
    {
        type            cyclic;

    }
    top
    {
	type zeroGradient;
    }
    bottom
    {
        type            fixedValue;
        value           uniform (0.01 0 0);
    }
    frontAndBack
    {
        type            empty;
    }
}


// ************************************************************************* //
