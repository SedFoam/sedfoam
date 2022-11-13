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
    location    "200";
    object      alpha.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar> 
200
(
0.400513
0.400581
0.400596
0.400658
0.400688
0.400741
0.400779
0.400828
0.400869
0.40092
0.40096
0.401013
0.401053
0.401109
0.401149
0.401207
0.401247
0.401306
0.401348
0.401408
0.401452
0.401513
0.401558
0.40162
0.401668
0.40173
0.40178
0.401843
0.401895
0.401959
0.402013
0.402079
0.402135
0.402202
0.40226
0.402329
0.402389
0.40246
0.402523
0.402595
0.40266
0.402735
0.402802
0.402879
0.40295
0.403029
0.403102
0.403184
0.40326
0.403346
0.403425
0.403513
0.403596
0.403688
0.403775
0.403871
0.403961
0.404061
0.404156
0.404261
0.40436
0.404471
0.404575
0.404691
0.404801
0.404924
0.40504
0.40517
0.405294
0.405431
0.405563
0.40571
0.405851
0.406009
0.40616
0.40633
0.406493
0.406677
0.406854
0.407055
0.40725
0.407471
0.407686
0.407932
0.408173
0.40845
0.408723
0.40904
0.409356
0.409727
0.410102
0.41055
0.411011
0.411579
0.412181
0.412955
0.413808
0.41502
0.41685
0.660811
0.999639
0.999999
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
        value           uniform 0.400513;
    }
    frontAndBackPlanes
    {
        type            empty;
    }
}


// ************************************************************************* //