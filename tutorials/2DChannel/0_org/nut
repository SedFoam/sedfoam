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
    object      nut;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0; 

boundaryField
{
    inletTop
    {
        type            calculated;
        value           uniform 0;
    }
    outletTop
    {
        type            calculated;
        value           uniform 0;
    }
    inlet
    {
        type            calculated;
        value           uniform 0;
    }
    outlet
    {
        type            calculated;
        value           uniform 0;
    }
    top
    {
//        type            symmetryPlane;
        type            calculated;
        value           uniform 0;
    }
    walls
    {
        type            calculated;
        value           uniform 0; 
       //type            zeroGradient;
    }
    apron
    {
        type            calculated;
        value           uniform 0; 
       //type            zeroGradient;
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
