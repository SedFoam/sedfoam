/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.1.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      epsilon;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

internalField   uniform 1e0; 

boundaryField
{
    inlet
    {
        type            fixedValue;
        value        	uniform  1;
    }
    outlet
    {
        type            zeroGradient;
    }
    top
    {
        type            zeroGradient;
    }
    walls
    {
        kn              0.0002;
        type            wilcoxOmegaWallFunction;
        value           uniform 594.884;
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
