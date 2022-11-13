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
    location    "1660";
    object      phi.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 3 -1 0 0 0 0];

oriented        oriented;

internalField   nonuniform List<scalar> 
119
(
5.41945e-17
1.09623e-16
1.64811e-16
2.21344e-16
2.77558e-16
3.3514e-16
3.92326e-16
4.50841e-16
5.08947e-16
5.68271e-16
6.27248e-16
6.87247e-16
7.47047e-16
8.07585e-16
8.68165e-16
9.29106e-16
9.90425e-16
1.05164e-15
1.11366e-15
1.17502e-15
1.23772e-15
1.29913e-15
1.36248e-15
1.42386e-15
1.48785e-15
1.54914e-15
1.61378e-15
1.67495e-15
1.74026e-15
1.80133e-15
1.86736e-15
1.92837e-15
1.99519e-15
2.05627e-15
2.12397e-15
2.18525e-15
2.25396e-15
2.31564e-15
2.38552e-15
2.44785e-15
2.51909e-15
2.58234e-15
2.65514e-15
2.71965e-15
2.79423e-15
2.86034e-15
2.93693e-15
3.00498e-15
3.08381e-15
3.15416e-15
3.23543e-15
3.3084e-15
3.39227e-15
3.46813e-15
3.55475e-15
3.63371e-15
3.72314e-15
3.80535e-15
3.89757e-15
3.98308e-15
4.07803e-15
4.1668e-15
4.26434e-15
4.3562e-15
4.45611e-15
4.55079e-15
4.65284e-15
4.74991e-15
4.85384e-15
4.95266e-15
5.0583e-15
5.15796e-15
5.26529e-15
5.36446e-15
5.47381e-15
5.57047e-15
5.68276e-15
5.77391e-15
5.89098e-15
5.97221e-15
6.09696e-15
6.16252e-15
6.29767e-15
6.3434e-15
6.48391e-15
6.52397e-15
6.6198e-15
6.78244e-15
6.50849e-15
1.11074e-14
5.12151e-15
2.09193e-18
2.8813e-17
-1.91036e-21
-2.41497e-21
-4.60574e-21
-5.00713e-22
6.39807e-21
7.21653e-21
3.51938e-21
1.1687e-21
-1.50051e-21
-5.25327e-21
-4.61259e-21
-1.04665e-21
-1.0477e-21
-8.32675e-22
4.23661e-21
7.01689e-21
2.42923e-21
-4.22853e-21
-5.22101e-21
-4.43801e-22
7.43278e-22
-1.44354e-21
3.5392e-22
3.37439e-21
3.6889e-21
-1.86364e-21
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
-4.151e-27
-1.24593e-26
-2.07826e-26
-2.91303e-26
-3.75127e-26
-4.59394e-26
-5.4421e-26
-6.2965e-26
-7.15825e-26
-8.02852e-26
-8.90799e-26
-9.79758e-26
-1.06986e-25
-1.16117e-25
-1.25378e-25
-1.34779e-25
-1.44328e-25
-1.54035e-25
-1.63908e-25
-1.73955e-25
-1.84185e-25
-1.94605e-25
-2.05223e-25
-2.16035e-25
-2.27082e-25
-2.38337e-25
-2.49822e-25
-2.61525e-25
-2.73469e-25
-2.85667e-25
-2.98095e-25
-3.10754e-25
-3.23677e-25
-3.36851e-25
-3.50276e-25
-3.63935e-25
-3.7753e-25
-3.91409e-25
-4.05463e-25
-4.19853e-25
-4.34366e-25
-4.49209e-25
-4.64217e-25
-4.79462e-25
-4.94892e-25
-5.10609e-25
-5.26396e-25
-5.42423e-25
-5.5851e-25
-5.7486e-25
-5.91284e-25
-6.07812e-25
-6.24425e-25
-6.411e-25
-6.57788e-25
-6.74548e-25
-6.91272e-25
-7.07958e-25
-7.24646e-25
-7.41328e-25
-7.57953e-25
-7.74429e-25
-7.90727e-25
-8.06816e-25
-8.22735e-25
-8.38168e-25
-8.53516e-25
-8.68457e-25
-8.83031e-25
-8.97207e-25
-9.11033e-25
-9.24236e-25
-9.37027e-25
-9.49296e-25
-9.61016e-25
-9.72157e-25
-9.82692e-25
-9.92597e-25
-1.00185e-24
-1.01041e-24
-1.01828e-24
-1.02543e-24
-1.03184e-24
-1.0375e-24
-1.04238e-24
-1.04669e-24
-1.04978e-24
-1.05228e-24
-1.05397e-24
-1.05483e-24
-1.05496e-24
-1.04322e-24
-1.0177e-24
-9.91588e-25
-9.64827e-25
-9.37362e-25
-9.09139e-25
-8.80109e-25
-8.50228e-25
-8.19457e-25
-7.87763e-25
-7.55121e-25
-7.21512e-25
-6.86926e-25
-6.51361e-25
-6.14821e-25
-5.7732e-25
-5.3888e-25
-4.99532e-25
-4.59314e-25
-4.18272e-25
-3.76459e-25
-3.33936e-25
-2.90771e-25
-2.47035e-25
-2.02808e-25
-1.58172e-25
-1.13215e-25
-6.80251e-26
-2.26952e-26
)
;
    }
    outlet
    {
        type            cyclic;
        value           nonuniform List<scalar> 
120
(
4.151e-27
1.24593e-26
2.07826e-26
2.91303e-26
3.75127e-26
4.59394e-26
5.4421e-26
6.2965e-26
7.15825e-26
8.02852e-26
8.90799e-26
9.79758e-26
1.06986e-25
1.16117e-25
1.25378e-25
1.34779e-25
1.44328e-25
1.54035e-25
1.63908e-25
1.73955e-25
1.84185e-25
1.94605e-25
2.05223e-25
2.16035e-25
2.27082e-25
2.38337e-25
2.49822e-25
2.61525e-25
2.73469e-25
2.85667e-25
2.98095e-25
3.10754e-25
3.23677e-25
3.36851e-25
3.50276e-25
3.63935e-25
3.7753e-25
3.91409e-25
4.05463e-25
4.19853e-25
4.34366e-25
4.49209e-25
4.64217e-25
4.79462e-25
4.94892e-25
5.10609e-25
5.26396e-25
5.42423e-25
5.5851e-25
5.7486e-25
5.91284e-25
6.07812e-25
6.24425e-25
6.411e-25
6.57788e-25
6.74548e-25
6.91272e-25
7.07958e-25
7.24646e-25
7.41328e-25
7.57953e-25
7.74429e-25
7.90727e-25
8.06816e-25
8.22735e-25
8.38168e-25
8.53516e-25
8.68457e-25
8.83031e-25
8.97207e-25
9.11033e-25
9.24236e-25
9.37027e-25
9.49296e-25
9.61016e-25
9.72157e-25
9.82692e-25
9.92597e-25
1.00185e-24
1.01041e-24
1.01828e-24
1.02543e-24
1.03184e-24
1.0375e-24
1.04238e-24
1.04669e-24
1.04978e-24
1.05228e-24
1.05397e-24
1.05483e-24
1.05496e-24
1.04322e-24
1.0177e-24
9.91588e-25
9.64827e-25
9.37362e-25
9.09139e-25
8.80109e-25
8.50228e-25
8.19457e-25
7.87763e-25
7.55121e-25
7.21512e-25
6.86926e-25
6.51361e-25
6.14821e-25
5.7732e-25
5.3888e-25
4.99532e-25
4.59314e-25
4.18272e-25
3.76459e-25
3.33936e-25
2.90771e-25
2.47035e-25
2.02808e-25
1.58172e-25
1.13215e-25
6.80251e-26
2.26952e-26
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