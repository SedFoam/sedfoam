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
    class       surfaceScalarField;
    location    "1760";
    object      phi.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 3 -1 0 0 0 0];

oriented        oriented;

internalField   nonuniform List<scalar> 
119
(
1.25628e-17
2.53571e-17
3.81996e-17
5.12243e-17
6.43787e-17
7.76485e-17
9.11265e-17
1.04628e-16
1.18448e-16
1.32165e-16
1.46348e-16
1.60261e-16
1.74832e-16
1.88919e-16
2.03902e-16
2.1814e-16
2.33561e-16
2.47928e-16
2.63812e-16
2.78285e-16
2.94657e-16
3.09212e-16
3.26096e-16
3.40714e-16
3.5813e-16
3.7279e-16
3.90758e-16
4.05444e-16
4.23979e-16
4.3868e-16
4.57789e-16
4.72498e-16
4.92187e-16
5.06903e-16
5.27169e-16
5.41894e-16
5.6273e-16
5.77474e-16
5.98864e-16
6.13642e-16
6.35563e-16
6.50394e-16
6.72816e-16
6.87728e-16
7.10615e-16
7.25635e-16
7.48946e-16
7.64108e-16
7.87794e-16
8.03132e-16
8.27143e-16
8.4269e-16
8.6697e-16
8.82758e-16
9.07253e-16
9.23303e-16
9.47961e-16
9.64285e-16
9.8906e-16
1.00565e-15
1.03051e-15
1.04732e-15
1.07226e-15
1.08922e-15
1.11424e-15
1.1312e-15
1.15641e-15
1.17313e-15
1.19867e-15
1.2148e-15
1.24094e-15
1.25594e-15
1.2831e-15
1.29625e-15
1.32499e-15
1.33541e-15
1.36632e-15
1.37318e-15
1.40644e-15
1.40975e-15
1.44368e-15
1.44683e-15
1.47346e-15
1.49082e-15
1.48247e-15
1.56393e-15
1.42894e-15
1.74819e-15
1.15263e-15
3.2803e-15
1.83448e-15
7.3635e-19
2.87256e-17
-3.7231e-21
-2.98918e-21
-1.41601e-21
1.04927e-21
3.12355e-21
3.79769e-21
3.20294e-21
1.26765e-21
-1.71905e-21
-3.6148e-21
-3.01741e-21
-1.4853e-21
-1.79666e-22
1.44612e-21
2.60987e-21
2.39959e-21
1.367e-21
-6.53701e-22
-2.81306e-21
-3.08159e-21
-1.52774e-21
1.14863e-21
3.11614e-21
2.98144e-21
8.22692e-22
-1.9088e-21
)
;

boundaryField
{
    inlet
    {
        type            cyclic;
        value           nonuniform List<scalar> 
120
(
-4.57224e-27
-1.37224e-26
-2.28859e-26
-3.20714e-26
-4.1288e-26
-5.05443e-26
-5.98492e-26
-6.92107e-26
-7.86378e-26
-8.81393e-26
-9.77227e-26
-1.07396e-25
-1.17168e-25
-1.27046e-25
-1.37037e-25
-1.47148e-25
-1.57387e-25
-1.67761e-25
-1.78275e-25
-1.88936e-25
-1.99749e-25
-2.1072e-25
-2.21854e-25
-2.33153e-25
-2.44628e-25
-2.56275e-25
-2.68101e-25
-2.80105e-25
-2.92292e-25
-3.04665e-25
-3.1722e-25
-3.29954e-25
-3.42875e-25
-3.55978e-25
-3.6926e-25
-3.82714e-25
-3.96011e-25
-4.09485e-25
-4.23104e-25
-4.36903e-25
-4.50823e-25
-4.6491e-25
-4.79116e-25
-4.93449e-25
-5.07891e-25
-5.22457e-25
-5.37084e-25
-5.51804e-25
-5.66562e-25
-5.81399e-25
-5.96255e-25
-6.11126e-25
-6.25996e-25
-6.40845e-25
-6.55648e-25
-6.70408e-25
-6.85084e-25
-6.99665e-25
-7.14147e-25
-7.28628e-25
-7.42967e-25
-7.57123e-25
-7.71071e-25
-7.84789e-25
-7.98269e-25
-8.11417e-25
-8.24319e-25
-8.36878e-25
-8.4909e-25
-8.60932e-25
-8.72402e-25
-8.83418e-25
-8.94017e-25
-9.04161e-25
-9.13827e-25
-9.22997e-25
-9.31652e-25
-9.39773e-25
-9.47345e-25
-9.5435e-25
-9.60774e-25
-9.66602e-25
-9.71822e-25
-9.76422e-25
-9.80392e-25
-9.83772e-25
-9.86404e-25
-9.88432e-25
-9.898e-25
-9.90504e-25
-9.90524e-25
-9.7692e-25
-9.47692e-25
-9.18043e-25
-8.87982e-25
-8.57506e-25
-8.26613e-25
-7.95303e-25
-7.63576e-25
-7.31433e-25
-6.98879e-25
-6.65917e-25
-6.32556e-25
-5.98802e-25
-5.64666e-25
-5.3016e-25
-4.95297e-25
-4.60093e-25
-4.24565e-25
-3.88731e-25
-3.52613e-25
-3.16232e-25
-2.79611e-25
-2.42776e-25
-2.05752e-25
-1.68567e-25
-1.31248e-25
-9.38243e-26
-5.63257e-26
-1.87817e-26
)
;
    }
    outlet
    {
        type            cyclic;
        value           nonuniform List<scalar> 
120
(
4.57224e-27
1.37224e-26
2.28859e-26
3.20714e-26
4.1288e-26
5.05443e-26
5.98492e-26
6.92107e-26
7.86378e-26
8.81393e-26
9.77227e-26
1.07396e-25
1.17168e-25
1.27046e-25
1.37037e-25
1.47148e-25
1.57387e-25
1.67761e-25
1.78275e-25
1.88936e-25
1.99749e-25
2.1072e-25
2.21854e-25
2.33153e-25
2.44628e-25
2.56275e-25
2.68101e-25
2.80105e-25
2.92292e-25
3.04665e-25
3.1722e-25
3.29954e-25
3.42875e-25
3.55978e-25
3.6926e-25
3.82714e-25
3.96011e-25
4.09485e-25
4.23104e-25
4.36903e-25
4.50823e-25
4.6491e-25
4.79116e-25
4.93449e-25
5.07891e-25
5.22457e-25
5.37084e-25
5.51804e-25
5.66562e-25
5.81399e-25
5.96255e-25
6.11126e-25
6.25996e-25
6.40845e-25
6.55648e-25
6.70408e-25
6.85084e-25
6.99665e-25
7.14147e-25
7.28628e-25
7.42967e-25
7.57123e-25
7.71071e-25
7.84789e-25
7.98269e-25
8.11417e-25
8.24319e-25
8.36878e-25
8.4909e-25
8.60932e-25
8.72402e-25
8.83418e-25
8.94017e-25
9.04161e-25
9.13827e-25
9.22997e-25
9.31652e-25
9.39773e-25
9.47345e-25
9.5435e-25
9.60774e-25
9.66602e-25
9.71822e-25
9.76422e-25
9.80392e-25
9.83772e-25
9.86404e-25
9.88432e-25
9.898e-25
9.90504e-25
9.90524e-25
9.7692e-25
9.47692e-25
9.18043e-25
8.87982e-25
8.57506e-25
8.26613e-25
7.95303e-25
7.63576e-25
7.31433e-25
6.98879e-25
6.65917e-25
6.32556e-25
5.98802e-25
5.64666e-25
5.3016e-25
4.95297e-25
4.60093e-25
4.24565e-25
3.88731e-25
3.52613e-25
3.16232e-25
2.79611e-25
2.42776e-25
2.05752e-25
1.68567e-25
1.31248e-25
9.38243e-26
5.63257e-26
1.87817e-26
)
;
    }
    top
    {
        type            fixedValue;
        value           uniform 0;
    }
    bottom
    {
        type            fixedValue;
        value           uniform 0;
    }
    frontAndBackPlanes
    {
        type            empty;
        value           nonuniform List<scalar> 0();
    }
}


// ************************************************************************* //