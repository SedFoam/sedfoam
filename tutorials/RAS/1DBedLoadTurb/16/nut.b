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
    location    "16";
    object      nut.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   nonuniform List<scalar> 
120
(
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
1.43641e-16
2.41e-15
1.68829e-14
8.6395e-14
4.02551e-13
3.35172e-12
2.85845e-11
1.45872e-10
4.91585e-10
1.31172e-09
3.26628e-09
8.44537e-09
2.35262e-08
6.94892e-08
2.08842e-07
6.19588e-07
1.82386e-06
5.31326e-06
1.42795e-05
3.30388e-05
6.51445e-05
0.000112107
0.000173409
0.000247342
0.00033183
0.000424882
0.000524726
0.000629791
0.000738627
0.000849848
0.000962103
0.00107412
0.00118484
0.00129359
0.00140022
0.00150515
0.00160906
0.00171263
0.00181618
0.00191969
0.00202286
0.00212527
0.00222651
0.00232619
0.00242402
0.00251975
0.00261321
0.00270423
0.0027927
0.00287851
0.00296159
0.00304186
0.00311924
0.00319365
0.00326505
0.00333334
0.00339846
0.00346036
0.00351895
0.00357415
0.00362588
0.00367405
0.00371856
0.00375932
0.00379622
0.00382915
0.00385799
0.0038826
0.00390284
0.00391857
0.00392962
0.00393581
0.00393693
0.00393279
0.00392313
0.0039077
0.00388621
0.00385833
0.00382369
0.00378188
0.00373243
0.00367478
0.00360832
0.00353228
0.0034458
0.00334778
0.00323692
0.00311156
0.00296958
0.00280811
0.00262329
0.00240908
0.00215816
0.00184432
0.00149416
0.000669548
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
        type            fixedValue;
        value           uniform 0;
    }
    frontAndBackPlanes
    {
        type            empty;
    }
}


// ************************************************************************* //