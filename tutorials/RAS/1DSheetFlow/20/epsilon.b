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
    object      epsilon.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -3 0 0 0 0];

internalField   nonuniform List<scalar> 
100
(
9.29851e-09
8.15186e-09
6.61882e-09
4.95163e-09
3.45579e-09
2.27638e-09
1.42502e-09
8.5068e-10
4.85607e-10
2.65493e-10
1.39268e-10
7.01646e-11
3.39913e-11
1.58445e-11
7.11027e-12
3.07277e-12
1.27857e-12
5.12242e-13
1.9742e-13
7.32747e-14
2.66613e-14
1.18321e-14
1.73126e-14
1.01354e-13
9.11151e-13
1.00718e-11
1.38073e-10
2.38537e-09
5.29523e-08
1.51906e-06
5.63099e-05
0.00184279
0.0105788
0.010605
0.00844742
0.0066708
0.00537143
0.00441742
0.00369902
0.00314323
0.00270239
0.00234575
0.00205135
0.00180501
0.00159615
0.00141658
0.00126109
0.00112488
0.00100463
0.000897764
0.000802212
0.000716411
0.000638906
0.000568581
0.000504514
0.000445941
0.000392226
0.00034284
0.000297741
0.000256547
0.000218519
0.000183428
0.000151117
0.000121973
9.59416e-05
7.28866e-05
5.2434e-05
3.47294e-05
2.00474e-05
9.13577e-06
2.36652e-06
2.42684e-07
2.68522e-08
1.01166e-08
7.69645e-09
7.23233e-09
7.1274e-09
7.09825e-09
7.08732e-09
7.08173e-09
7.07889e-09
7.07867e-09
7.08162e-09
7.08868e-09
7.09676e-09
7.09777e-09
7.09167e-09
7.07786e-09
7.07082e-09
7.0705e-09
7.05954e-09
7.05458e-09
7.04837e-09
7.04777e-09
7.05339e-09
7.0666e-09
7.0662e-09
7.07633e-09
7.07366e-09
7.08444e-09
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
        type            zeroGradient;
    }
    bottom
    {
        type            zeroGradient;
    }
    frontAndBackPlanes
    {
        type            empty;
    }
}


// ************************************************************************* //