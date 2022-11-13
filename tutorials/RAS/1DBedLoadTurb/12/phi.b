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
    location    "12";
    object      phi.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 3 -1 0 0 0 0];

oriented        oriented;

internalField   nonuniform List<scalar> 
119
(
-1.68311e-14
-4.37204e-14
-6.26059e-14
-8.87823e-14
-1.11043e-13
-1.3837e-13
-1.65395e-13
-1.94475e-13
-2.25272e-13
-2.58148e-13
-2.93256e-13
-3.30991e-13
-3.71636e-13
-4.15847e-13
-4.64173e-13
-5.17385e-13
-5.7681e-13
-6.44103e-13
-7.2139e-13
-8.10545e-13
-9.12671e-13
-1.02865e-12
-1.16107e-12
-1.32093e-12
-1.53188e-12
-1.82162e-12
-2.25941e-12
-3.08045e-12
-4.73497e-12
-5.43826e-12
-6.31743e-12
-7.58222e-12
-9.64894e-12
-1.26961e-11
-1.76334e-11
-2.64261e-11
-4.11722e-11
-6.27684e-11
-9.0071e-11
-1.21065e-10
-1.53446e-10
-1.84524e-10
-2.11749e-10
-2.33413e-10
-2.48846e-10
-2.58184e-10
-2.62032e-10
-2.6117e-10
-2.56347e-10
-2.48171e-10
-2.37076e-10
-2.23335e-10
-2.07106e-10
-1.88515e-10
-1.67756e-10
-1.45233e-10
-1.2169e-10
-9.82567e-11
-7.62993e-11
-5.70592e-11
-4.12809e-11
-2.90772e-11
-2.00729e-11
-1.3661e-11
-9.21151e-12
-6.1802e-12
-4.14094e-12
-2.77807e-12
-1.86989e-12
-1.24043e-12
-8.38412e-13
-4.95116e-13
-3.17188e-13
-1.46946e-13
-4.61894e-14
4.60338e-15
3.64965e-14
7.65378e-14
2.43037e-15
1.96343e-12
-1.22762e-16
-1.30881e-16
-1.0169e-16
-4.40755e-17
-1.11536e-16
-5.08021e-17
-1.56068e-16
-1.28665e-16
-9.14435e-17
-1.17884e-16
-1.76873e-16
-1.92989e-16
-1.56536e-16
-1.4014e-16
-1.5579e-16
-1.86733e-16
-2.29788e-16
-2.25187e-16
-2.25964e-16
-2.19437e-16
-2.24527e-16
-2.49684e-16
-2.57989e-16
-2.84852e-16
-2.67894e-16
-2.87851e-16
-2.92832e-16
-3.0355e-16
-3.20074e-16
-3.25837e-16
-3.28486e-16
-3.37774e-16
-3.45877e-16
-3.60933e-16
-3.74265e-16
-3.72768e-16
-3.71544e-16
-3.77528e-16
-3.81202e-16
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
-3.77788e-08
-3.83823e-08
-3.84032e-08
-3.84318e-08
-3.8455e-08
-3.84813e-08
-3.85073e-08
-3.85342e-08
-3.85616e-08
-3.85898e-08
-3.86187e-08
-3.86484e-08
-3.86789e-08
-3.87103e-08
-3.87427e-08
-3.87761e-08
-3.88106e-08
-3.88463e-08
-3.88834e-08
-3.8922e-08
-3.89624e-08
-3.90052e-08
-3.90508e-08
-3.9099e-08
-3.91488e-08
-3.92003e-08
-3.92561e-08
-3.93262e-08
-3.94416e-08
-3.96249e-08
-3.99925e-08
-4.18377e-08
-4.76413e-08
-5.82308e-08
-7.28309e-08
-9.12518e-08
-1.14277e-07
-1.4329e-07
-1.80429e-07
-2.28839e-07
-2.93397e-07
-3.83965e-07
-5.19625e-07
-7.16901e-07
-9.73182e-07
-1.27097e-06
-1.58924e-06
-1.91077e-06
-2.22416e-06
-2.52298e-06
-2.80418e-06
-3.0667e-06
-3.31056e-06
-3.53627e-06
-3.74454e-06
-3.93616e-06
-4.11194e-06
-4.27278e-06
-4.41969e-06
-4.55382e-06
-4.67641e-06
-4.78875e-06
-4.89205e-06
-4.98742e-06
-5.07582e-06
-5.15806e-06
-5.23483e-06
-5.30674e-06
-5.37427e-06
-5.43786e-06
-5.49786e-06
-5.55459e-06
-5.60834e-06
-5.65933e-06
-5.70779e-06
-5.75389e-06
-5.79782e-06
-5.83971e-06
-5.8797e-06
-5.9179e-06
-5.95444e-06
-5.98939e-06
-6.02287e-06
-6.05493e-06
-6.08567e-06
-6.11514e-06
-6.14341e-06
-6.17054e-06
-6.19657e-06
-6.22156e-06
-6.24555e-06
-6.26858e-06
-6.29069e-06
-6.31191e-06
-6.33228e-06
-6.35182e-06
-6.37056e-06
-6.38853e-06
-6.40575e-06
-6.42224e-06
-6.43802e-06
-6.45311e-06
-6.46752e-06
-6.48128e-06
-6.49438e-06
-6.50684e-06
-6.51867e-06
-6.52988e-06
-6.54047e-06
-6.55044e-06
-6.5598e-06
-6.56853e-06
-6.57664e-06
-6.58411e-06
-6.59093e-06
-6.59706e-06
-6.60246e-06
-6.60709e-06
-6.6108e-06
-6.61366e-06
)
;
    }
    outlet
    {
        type            cyclic;
        value           nonuniform List<scalar> 
120
(
3.77788e-08
3.83823e-08
3.84032e-08
3.84318e-08
3.8455e-08
3.84813e-08
3.85073e-08
3.85342e-08
3.85616e-08
3.85898e-08
3.86187e-08
3.86484e-08
3.86789e-08
3.87103e-08
3.87427e-08
3.87761e-08
3.88106e-08
3.88463e-08
3.88834e-08
3.8922e-08
3.89624e-08
3.90052e-08
3.90508e-08
3.9099e-08
3.91488e-08
3.92003e-08
3.92561e-08
3.93262e-08
3.94416e-08
3.96249e-08
3.99925e-08
4.18377e-08
4.76413e-08
5.82308e-08
7.28309e-08
9.12518e-08
1.14277e-07
1.4329e-07
1.80429e-07
2.28839e-07
2.93397e-07
3.83965e-07
5.19625e-07
7.16901e-07
9.73182e-07
1.27097e-06
1.58924e-06
1.91077e-06
2.22416e-06
2.52298e-06
2.80418e-06
3.0667e-06
3.31056e-06
3.53627e-06
3.74454e-06
3.93616e-06
4.11194e-06
4.27278e-06
4.41969e-06
4.55382e-06
4.67641e-06
4.78875e-06
4.89205e-06
4.98742e-06
5.07582e-06
5.15806e-06
5.23483e-06
5.30674e-06
5.37427e-06
5.43786e-06
5.49786e-06
5.55459e-06
5.60834e-06
5.65933e-06
5.70779e-06
5.75389e-06
5.79782e-06
5.83971e-06
5.8797e-06
5.9179e-06
5.95444e-06
5.98939e-06
6.02287e-06
6.05493e-06
6.08567e-06
6.11514e-06
6.14341e-06
6.17054e-06
6.19657e-06
6.22156e-06
6.24555e-06
6.26858e-06
6.29069e-06
6.31191e-06
6.33228e-06
6.35182e-06
6.37056e-06
6.38853e-06
6.40575e-06
6.42224e-06
6.43802e-06
6.45311e-06
6.46752e-06
6.48128e-06
6.49438e-06
6.50684e-06
6.51867e-06
6.52988e-06
6.54047e-06
6.55044e-06
6.5598e-06
6.56853e-06
6.57664e-06
6.58411e-06
6.59093e-06
6.59706e-06
6.60246e-06
6.60709e-06
6.6108e-06
6.61366e-06
)
;
    }
    top
    {
        type            calculated;
        value           uniform -3.94893e-16;
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