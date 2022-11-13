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
    location    "14";
    object      phi.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 3 -1 0 0 0 0];

oriented        oriented;

internalField   nonuniform List<scalar> 
119
(
-1.13958e-14
-3.10477e-14
-4.30356e-14
-6.21715e-14
-7.69875e-14
-9.63434e-14
-1.14594e-13
-1.34571e-13
-1.55679e-13
-1.78265e-13
-2.02409e-13
-2.28381e-13
-2.56555e-13
-2.87309e-13
-3.21021e-13
-3.58396e-13
-4.00558e-13
-4.48824e-13
-5.04784e-13
-5.70397e-13
-6.47333e-13
-7.36719e-13
-8.41144e-13
-9.69484e-13
-1.14012e-12
-1.38182e-12
-1.77095e-12
-2.44612e-12
-3.11331e-12
-3.46786e-12
-3.97756e-12
-4.76652e-12
-6.02823e-12
-7.89718e-12
-1.10435e-11
-1.66368e-11
-2.56085e-11
-3.82478e-11
-5.36958e-11
-7.08364e-11
-8.8412e-11
-1.04944e-10
-1.19144e-10
-1.30221e-10
-1.37919e-10
-1.42373e-10
-1.43936e-10
-1.43028e-10
-1.40047e-10
-1.3531e-10
-1.29043e-10
-1.21391e-10
-1.12436e-10
-1.02249e-10
-9.09333e-11
-7.87067e-11
-6.597e-11
-5.33443e-11
-4.16025e-11
-3.14688e-11
-2.3384e-11
-1.73957e-11
-1.3219e-11
-1.03892e-11
-8.42555e-12
-6.91625e-12
-5.56605e-12
-4.21319e-12
-2.82559e-12
-1.47549e-12
-3.01472e-13
5.6671e-13
1.07672e-12
1.26776e-12
1.21999e-12
1.15338e-12
9.50452e-13
1.04023e-12
6.69238e-13
2.40577e-12
-4.73388e-17
-8.67335e-17
-1.13018e-16
-1.25997e-16
-1.42415e-16
-1.44233e-16
-1.24324e-16
-8.85032e-17
-2.761e-17
9.48482e-18
2.59637e-17
-2.84288e-17
1.23376e-16
9.74113e-17
1.17717e-16
1.58899e-16
1.32537e-16
8.44003e-17
7.86831e-17
9.56043e-17
7.79106e-17
3.73623e-17
2.30445e-17
1.82693e-17
1.73688e-17
-6.34546e-18
-6.51006e-17
-5.99262e-17
-5.56083e-17
-6.95871e-17
-6.43668e-17
-6.59192e-17
-6.75947e-17
-5.18896e-17
-5.29654e-17
-8.76512e-17
-1.11874e-16
-1.12015e-16
-9.68078e-17
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
-3.77795e-08
-3.83834e-08
-3.84039e-08
-3.84329e-08
-3.84559e-08
-3.84824e-08
-3.85085e-08
-3.85354e-08
-3.8563e-08
-3.85913e-08
-3.86203e-08
-3.86501e-08
-3.86807e-08
-3.87123e-08
-3.87449e-08
-3.87785e-08
-3.88133e-08
-3.88493e-08
-3.88868e-08
-3.89259e-08
-3.8967e-08
-3.90104e-08
-3.90568e-08
-3.91062e-08
-3.91582e-08
-3.9213e-08
-3.92756e-08
-3.93618e-08
-3.95035e-08
-3.96917e-08
-4.04019e-08
-4.35077e-08
-5.11342e-08
-6.32582e-08
-7.91409e-08
-9.89408e-08
-1.23637e-07
-1.54768e-07
-1.94687e-07
-2.46792e-07
-3.16517e-07
-4.15196e-07
-5.62307e-07
-7.70767e-07
-1.03457e-06
-1.33554e-06
-1.6537e-06
-1.9733e-06
-2.28407e-06
-2.5803e-06
-2.8593e-06
-3.12017e-06
-3.36294e-06
-3.58812e-06
-3.79637e-06
-3.9884e-06
-4.16495e-06
-4.32683e-06
-4.47496e-06
-4.61038e-06
-4.73427e-06
-4.84785e-06
-4.95234e-06
-5.04881e-06
-5.13823e-06
-5.22142e-06
-5.2991e-06
-5.37184e-06
-5.44017e-06
-5.50449e-06
-5.5652e-06
-5.62261e-06
-5.67698e-06
-5.72858e-06
-5.77762e-06
-5.82427e-06
-5.86872e-06
-5.91111e-06
-5.95158e-06
-5.99025e-06
-6.02722e-06
-6.0626e-06
-6.09648e-06
-6.12893e-06
-6.16004e-06
-6.18987e-06
-6.21848e-06
-6.24594e-06
-6.27229e-06
-6.29758e-06
-6.32186e-06
-6.34517e-06
-6.36755e-06
-6.38903e-06
-6.40964e-06
-6.42942e-06
-6.44839e-06
-6.46657e-06
-6.484e-06
-6.50069e-06
-6.51666e-06
-6.53193e-06
-6.54651e-06
-6.56043e-06
-6.57369e-06
-6.58629e-06
-6.59826e-06
-6.6096e-06
-6.62031e-06
-6.6304e-06
-6.63986e-06
-6.6487e-06
-6.6569e-06
-6.66445e-06
-6.67134e-06
-6.67753e-06
-6.683e-06
-6.68768e-06
-6.69142e-06
-6.6943e-06
)
;
    }
    outlet
    {
        type            cyclic;
        value           nonuniform List<scalar> 
120
(
3.77795e-08
3.83834e-08
3.84039e-08
3.84329e-08
3.84559e-08
3.84824e-08
3.85085e-08
3.85354e-08
3.8563e-08
3.85913e-08
3.86203e-08
3.86501e-08
3.86807e-08
3.87123e-08
3.87449e-08
3.87785e-08
3.88133e-08
3.88493e-08
3.88868e-08
3.89259e-08
3.8967e-08
3.90104e-08
3.90568e-08
3.91062e-08
3.91582e-08
3.9213e-08
3.92756e-08
3.93618e-08
3.95035e-08
3.96917e-08
4.04019e-08
4.35077e-08
5.11342e-08
6.32582e-08
7.91409e-08
9.89408e-08
1.23637e-07
1.54768e-07
1.94687e-07
2.46792e-07
3.16517e-07
4.15196e-07
5.62307e-07
7.70767e-07
1.03457e-06
1.33554e-06
1.6537e-06
1.9733e-06
2.28407e-06
2.5803e-06
2.8593e-06
3.12017e-06
3.36294e-06
3.58812e-06
3.79637e-06
3.9884e-06
4.16495e-06
4.32683e-06
4.47496e-06
4.61038e-06
4.73427e-06
4.84785e-06
4.95234e-06
5.04881e-06
5.13823e-06
5.22142e-06
5.2991e-06
5.37184e-06
5.44017e-06
5.50449e-06
5.5652e-06
5.62261e-06
5.67698e-06
5.72858e-06
5.77762e-06
5.82427e-06
5.86872e-06
5.91111e-06
5.95158e-06
5.99025e-06
6.02722e-06
6.0626e-06
6.09648e-06
6.12893e-06
6.16004e-06
6.18987e-06
6.21848e-06
6.24594e-06
6.27229e-06
6.29758e-06
6.32186e-06
6.34517e-06
6.36755e-06
6.38903e-06
6.40964e-06
6.42942e-06
6.44839e-06
6.46657e-06
6.484e-06
6.50069e-06
6.51666e-06
6.53193e-06
6.54651e-06
6.56043e-06
6.57369e-06
6.58629e-06
6.59826e-06
6.6096e-06
6.62031e-06
6.6304e-06
6.63986e-06
6.6487e-06
6.6569e-06
6.66445e-06
6.67134e-06
6.67753e-06
6.683e-06
6.68768e-06
6.69142e-06
6.6943e-06
)
;
    }
    top
    {
        type            calculated;
        value           uniform -8.64059e-17;
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