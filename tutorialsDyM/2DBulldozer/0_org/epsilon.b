/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1906                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      delta;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -3 0 0 0 0];
internalField   uniform 1e-8; 

boundaryField
{
    hole
    {
        type            zeroGradient;
    }
    "overset.*"
    {
       type            overset;
    }
    inandouthalf21
    {
	type		zeroGradient;
    }
    inandouthalf22
    {
	type		zeroGradient;
    }
    top
    {
        type            zeroGradient;

    }
    bottom
    {
        type            zeroGradient;
    }
    inlet
    {
        type            zeroGradient;
    }
    outlet
    {
        type            zeroGradient;
    }
    walls
    {
        type            zeroGradient;
    }
    frontAndBackPlanes
    {
        type            empty;
    }
}


// ************************************************************************* //
