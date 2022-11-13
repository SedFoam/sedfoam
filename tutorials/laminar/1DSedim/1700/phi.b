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
    location    "1700";
    object      phi.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 3 -1 0 0 0 0];

oriented        oriented;

internalField   nonuniform List<scalar> 
119
(
3.02025e-17
6.10442e-17
9.18176e-17
1.2332e-16
1.54741e-16
1.8689e-16
2.18973e-16
2.51739e-16
2.84503e-16
3.17853e-16
3.51323e-16
3.85224e-16
4.1943e-16
4.53844e-16
4.8882e-16
5.23713e-16
5.59497e-16
5.94839e-16
6.31472e-16
6.67234e-16
7.04758e-16
7.4092e-16
7.7938e-16
8.15927e-16
8.55366e-16
8.9229e-16
9.32748e-16
9.70052e-16
1.01157e-15
1.04926e-15
1.09185e-15
1.12995e-15
1.17365e-15
1.21217e-15
1.25699e-15
1.29596e-15
1.34189e-15
1.38133e-15
1.42835e-15
1.4683e-15
1.51638e-15
1.55686e-15
1.60593e-15
1.64695e-15
1.69695e-15
1.73853e-15
1.78935e-15
1.83149e-15
1.88303e-15
1.9257e-15
1.97783e-15
2.02102e-15
2.07358e-15
2.11725e-15
2.17009e-15
2.21419e-15
2.26714e-15
2.31159e-15
2.36451e-15
2.40923e-15
2.46194e-15
2.50682e-15
2.5592e-15
2.60409e-15
2.65605e-15
2.70077e-15
2.75227e-15
2.79654e-15
2.84764e-15
2.89108e-15
2.94198e-15
2.98399e-15
3.03515e-15
3.07479e-15
3.12699e-15
3.16287e-15
3.21736e-15
3.24746e-15
3.3059e-15
3.32771e-15
3.39163e-15
3.40325e-15
3.47121e-15
3.47652e-15
3.533e-15
3.56399e-15
3.53075e-15
3.76491e-15
3.2087e-15
6.91135e-15
3.33066e-15
1.29945e-18
2.87597e-17
-6.85965e-21
-5.59907e-21
-2.88342e-21
1.58318e-21
5.36969e-21
6.77381e-21
6.25216e-21
3.14827e-21
-2.60587e-21
-6.47902e-21
-5.31358e-21
-2.70821e-21
-1.01657e-21
1.77068e-21
4.12758e-21
3.9713e-21
3.04442e-21
2.24844e-22
-4.06325e-21
-5.12105e-21
-3.27839e-21
7.23374e-22
4.76588e-21
4.78257e-21
1.83357e-21
-2.84972e-21
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
-4.32465e-27
-1.29801e-26
-2.16501e-26
-3.03438e-26
-3.90711e-26
-4.78413e-26
-5.66644e-26
-6.55486e-26
-7.45037e-26
-8.35404e-26
-9.26661e-26
-1.0189e-25
-1.11222e-25
-1.2067e-25
-1.30243e-25
-1.39949e-25
-1.49797e-25
-1.59794e-25
-1.69948e-25
-1.80267e-25
-1.90758e-25
-2.01427e-25
-2.12281e-25
-2.2332e-25
-2.34569e-25
-2.46013e-25
-2.57667e-25
-2.69525e-25
-2.81601e-25
-2.93902e-25
-3.06415e-25
-3.19141e-25
-3.32097e-25
-3.45276e-25
-3.58677e-25
-3.72287e-25
-3.85793e-25
-3.99535e-25
-4.13444e-25
-4.27616e-25
-4.41918e-25
-4.56473e-25
-4.71176e-25
-4.86068e-25
-5.01111e-25
-5.16364e-25
-5.31689e-25
-5.47185e-25
-5.62735e-25
-5.78458e-25
-5.94233e-25
-6.1007e-25
-6.25952e-25
-6.41859e-25
-6.57752e-25
-6.73662e-25
-6.89514e-25
-7.053e-25
-7.21039e-25
-7.36768e-25
-7.52398e-25
-7.67861e-25
-7.83131e-25
-7.98181e-25
-8.13025e-25
-8.27465e-25
-8.4173e-25
-8.55621e-25
-8.69152e-25
-8.82296e-25
-8.95073e-25
-9.07313e-25
-9.19132e-25
-9.30459e-25
-9.41267e-25
-9.51533e-25
-9.61233e-25
-9.70345e-25
-9.78849e-25
-9.86723e-25
-9.9395e-25
-1.00051e-24
-1.00639e-24
-1.01158e-24
-1.01606e-24
-1.01993e-24
-1.02284e-24
-1.02513e-24
-1.02668e-24
-1.02747e-24
-1.02754e-24
-1.01457e-24
-9.86574e-25
-9.58067e-25
-9.29041e-25
-8.99471e-25
-8.69335e-25
-8.38613e-25
-8.07288e-25
-7.75347e-25
-7.42779e-25
-7.09576e-25
-6.75738e-25
-6.41264e-25
-6.0616e-25
-5.70437e-25
-5.34107e-25
-4.97191e-25
-4.59711e-25
-4.21693e-25
-3.8317e-25
-3.44177e-25
-3.04751e-25
-2.64936e-25
-2.24777e-25
-1.84322e-25
-1.43621e-25
-1.02727e-25
-6.16938e-26
-2.05767e-26
)
;
    }
    outlet
    {
        type            cyclic;
        value           nonuniform List<scalar> 
120
(
4.32465e-27
1.29801e-26
2.16501e-26
3.03438e-26
3.90711e-26
4.78413e-26
5.66644e-26
6.55486e-26
7.45037e-26
8.35404e-26
9.26661e-26
1.0189e-25
1.11222e-25
1.2067e-25
1.30243e-25
1.39949e-25
1.49797e-25
1.59794e-25
1.69948e-25
1.80267e-25
1.90758e-25
2.01427e-25
2.12281e-25
2.2332e-25
2.34569e-25
2.46013e-25
2.57667e-25
2.69525e-25
2.81601e-25
2.93902e-25
3.06415e-25
3.19141e-25
3.32097e-25
3.45276e-25
3.58677e-25
3.72287e-25
3.85793e-25
3.99535e-25
4.13444e-25
4.27616e-25
4.41918e-25
4.56473e-25
4.71176e-25
4.86068e-25
5.01111e-25
5.16364e-25
5.31689e-25
5.47185e-25
5.62735e-25
5.78458e-25
5.94233e-25
6.1007e-25
6.25952e-25
6.41859e-25
6.57752e-25
6.73662e-25
6.89514e-25
7.053e-25
7.21039e-25
7.36768e-25
7.52398e-25
7.67861e-25
7.83131e-25
7.98181e-25
8.13025e-25
8.27465e-25
8.4173e-25
8.55621e-25
8.69152e-25
8.82296e-25
8.95073e-25
9.07313e-25
9.19132e-25
9.30459e-25
9.41267e-25
9.51533e-25
9.61233e-25
9.70345e-25
9.78849e-25
9.86723e-25
9.9395e-25
1.00051e-24
1.00639e-24
1.01158e-24
1.01606e-24
1.01993e-24
1.02284e-24
1.02513e-24
1.02668e-24
1.02747e-24
1.02754e-24
1.01457e-24
9.86574e-25
9.58067e-25
9.29041e-25
8.99471e-25
8.69335e-25
8.38613e-25
8.07288e-25
7.75347e-25
7.42779e-25
7.09576e-25
6.75738e-25
6.41264e-25
6.0616e-25
5.70437e-25
5.34107e-25
4.97191e-25
4.59711e-25
4.21693e-25
3.8317e-25
3.44177e-25
3.04751e-25
2.64936e-25
2.24777e-25
1.84322e-25
1.43621e-25
1.02727e-25
6.16938e-26
2.05767e-26
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