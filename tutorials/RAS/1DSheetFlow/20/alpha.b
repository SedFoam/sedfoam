/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2206                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       volScalarField;
    location    "20";
    object      alpha.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar> 
100
(
0.449325
0.449464
0.449455
0.449614
0.449616
0.44978
0.449795
0.449965
0.449992
0.450171
0.450212
0.450402
0.450459
0.450665
0.450739
0.450966
0.451059
0.451316
0.451432
0.451728
0.451874
0.452221
0.452419
0.452823
0.453127
0.453566
0.454164
0.454483
0.456162
0.456227
0.460947
0.642259
0.860234
0.949929
0.979073
0.99105
0.995873
0.998247
0.999133
0.999658
0.999811
0.999939
0.99996
0.999994
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
)
;

boundaryField
{
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
        type            calculated;
        value           uniform 1;
    }
    bottom
    {
        type            calculated;
        value           uniform 0.449325;
    }
    frontAndBackPlanes
    {
        type            empty;
    }
}


// ************************************************************************* //
