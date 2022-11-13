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
    location    "1280";
    object      alpha.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar> 
120
(
0.400514
0.400613
0.400612
0.400721
0.400724
0.400833
0.400838
0.400948
0.400955
0.401066
0.401076
0.401187
0.4012
0.401311
0.401328
0.401438
0.401459
0.401569
0.401594
0.401704
0.401733
0.401843
0.401876
0.401986
0.402024
0.402133
0.402177
0.402286
0.402335
0.402444
0.402499
0.402608
0.402668
0.402777
0.402844
0.402954
0.403026
0.403137
0.403216
0.403328
0.403413
0.403527
0.403619
0.403736
0.403834
0.403954
0.40406
0.404184
0.404297
0.404425
0.404546
0.404679
0.404808
0.404949
0.405087
0.405235
0.405382
0.405539
0.405698
0.405865
0.406036
0.406215
0.4064
0.406594
0.406795
0.407006
0.407227
0.407459
0.407703
0.407961
0.408235
0.408525
0.408836
0.409169
0.409528
0.409919
0.410349
0.410826
0.411367
0.412006
0.412788
0.413693
0.414755
0.416249
0.41761
0.419878
0.455983
0.477952
0.489235
0.501097
0.476312
0.73151
0.999533
0.999995
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
        value           uniform 0.400514;
    }
    frontAndBackPlanes
    {
        type            empty;
    }
}


// ************************************************************************* //