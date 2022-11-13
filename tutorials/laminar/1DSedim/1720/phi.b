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
    location    "1720";
    object      phi.b;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 3 -1 0 0 0 0];

oriented        oriented;

internalField   nonuniform List<scalar> 
119
(
2.36468e-17
4.77868e-17
7.18894e-17
9.6529e-17
1.21159e-16
1.46292e-16
1.71461e-16
1.97061e-16
2.22784e-16
2.48822e-16
2.75114e-16
3.01557e-16
3.2844e-16
3.55249e-16
3.8274e-16
4.09875e-16
4.37996e-16
4.65415e-16
4.94185e-16
5.21843e-16
5.51277e-16
5.7913e-16
6.09242e-16
6.3725e-16
6.68047e-16
6.96167e-16
7.27654e-16
7.55851e-16
7.88024e-16
8.16267e-16
8.49116e-16
8.77379e-16
9.10888e-16
9.39153e-16
9.73298e-16
1.00156e-15
1.03631e-15
1.06456e-15
1.09988e-15
1.12814e-15
1.16399e-15
1.19227e-15
1.22862e-15
1.25695e-15
1.29374e-15
1.32217e-15
1.35937e-15
1.38795e-15
1.42551e-15
1.4543e-15
1.4922e-15
1.52127e-15
1.55945e-15
1.58887e-15
1.62731e-15
1.65715e-15
1.69583e-15
1.72615e-15
1.76504e-15
1.79588e-15
1.83496e-15
1.86633e-15
1.90559e-15
1.93744e-15
1.9769e-15
2.0091e-15
2.04881e-15
2.08112e-15
2.12117e-15
2.15319e-15
2.19381e-15
2.22491e-15
2.26649e-15
2.29569e-15
2.33888e-15
2.36481e-15
2.41055e-15
2.43143e-15
2.48081e-15
2.49477e-15
2.54812e-15
2.5549e-15
2.60836e-15
2.61561e-15
2.64857e-15
2.6966e-15
2.62111e-15
2.89862e-15
2.29326e-15
5.42307e-15
2.72871e-15
1.07009e-18
2.87461e-17
-5.26174e-21
-4.36776e-21
-2.35988e-21
1.04248e-21
4.03897e-21
5.22848e-21
4.81572e-21
2.38789e-21
-1.93922e-21
-4.82969e-21
-4.03414e-21
-2.14431e-21
-8.73381e-22
1.16508e-21
2.96289e-21
3.0588e-21
2.44302e-21
2.82302e-22
-2.94463e-21
-3.92319e-21
-2.57261e-21
5.74318e-22
3.50055e-21
3.58926e-21
1.29672e-21
-2.03662e-21
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
-4.40921e-27
-1.32337e-26
-2.20723e-26
-3.09342e-26
-3.98289e-26
-4.87656e-26
-5.7754e-26
-6.68022e-26
-7.59197e-26
-8.51167e-26
-9.44005e-26
-1.0378e-25
-1.13265e-25
-1.22863e-25
-1.32582e-25
-1.42431e-25
-1.52416e-25
-1.62547e-25
-1.72829e-25
-1.83271e-25
-1.93878e-25
-2.04657e-25
-2.15614e-25
-2.26749e-25
-2.38083e-25
-2.49604e-25
-2.61326e-25
-2.73243e-25
-2.85366e-25
-2.97702e-25
-3.10239e-25
-3.22977e-25
-3.35932e-25
-3.49095e-25
-3.62466e-25
-3.76033e-25
-3.89478e-25
-4.03137e-25
-4.16957e-25
-4.31009e-25
-4.45191e-25
-4.59593e-25
-4.74134e-25
-4.88842e-25
-5.03686e-25
-5.18708e-25
-5.338e-25
-5.49033e-25
-5.64316e-25
-5.79736e-25
-5.95196e-25
-6.10702e-25
-6.26236e-25
-6.41777e-25
-6.57294e-25
-6.72803e-25
-6.88244e-25
-7.03607e-25
-7.18904e-25
-7.34197e-25
-7.49375e-25
-7.64379e-25
-7.79182e-25
-7.93761e-25
-8.08119e-25
-8.22106e-25
-8.35885e-25
-8.49302e-25
-8.62364e-25
-8.75045e-25
-8.87353e-25
-8.99159e-25
-9.10543e-25
-9.21447e-25
-9.31848e-25
-9.41722e-25
-9.51048e-25
-9.59807e-25
-9.67977e-25
-9.7554e-25
-9.8248e-25
-9.88779e-25
-9.94424e-25
-9.994e-25
-1.0037e-24
-1.00739e-24
-1.0102e-24
-1.0124e-24
-1.01388e-24
-1.01464e-24
-1.01469e-24
-1.0014e-24
-9.72768e-25
-9.43658e-25
-9.14073e-25
-8.83997e-25
-8.53418e-25
-8.22325e-25
-7.90708e-25
-7.58562e-25
-7.25881e-25
-6.92664e-25
-6.58914e-25
-6.24634e-25
-5.89833e-25
-5.54521e-25
-5.18712e-25
-4.82424e-25
-4.45678e-25
-4.08496e-25
-3.70906e-25
-3.32937e-25
-2.94621e-25
-2.55992e-25
-2.17088e-25
-1.77947e-25
-1.38611e-25
-9.91195e-26
-5.95176e-26
-1.98488e-26
)
;
    }
    outlet
    {
        type            cyclic;
        value           nonuniform List<scalar> 
120
(
4.40921e-27
1.32337e-26
2.20723e-26
3.09342e-26
3.98289e-26
4.87656e-26
5.7754e-26
6.68022e-26
7.59197e-26
8.51167e-26
9.44005e-26
1.0378e-25
1.13265e-25
1.22863e-25
1.32582e-25
1.42431e-25
1.52416e-25
1.62547e-25
1.72829e-25
1.83271e-25
1.93878e-25
2.04657e-25
2.15614e-25
2.26749e-25
2.38083e-25
2.49604e-25
2.61326e-25
2.73243e-25
2.85366e-25
2.97702e-25
3.10239e-25
3.22977e-25
3.35932e-25
3.49095e-25
3.62466e-25
3.76033e-25
3.89478e-25
4.03137e-25
4.16957e-25
4.31009e-25
4.45191e-25
4.59593e-25
4.74134e-25
4.88842e-25
5.03686e-25
5.18708e-25
5.338e-25
5.49033e-25
5.64316e-25
5.79736e-25
5.95196e-25
6.10702e-25
6.26236e-25
6.41777e-25
6.57294e-25
6.72803e-25
6.88244e-25
7.03607e-25
7.18904e-25
7.34197e-25
7.49375e-25
7.64379e-25
7.79182e-25
7.93761e-25
8.08119e-25
8.22106e-25
8.35885e-25
8.49302e-25
8.62364e-25
8.75045e-25
8.87353e-25
8.99159e-25
9.10543e-25
9.21447e-25
9.31848e-25
9.41722e-25
9.51048e-25
9.59807e-25
9.67977e-25
9.7554e-25
9.8248e-25
9.88779e-25
9.94424e-25
9.994e-25
1.0037e-24
1.00739e-24
1.0102e-24
1.0124e-24
1.01388e-24
1.01464e-24
1.01469e-24
1.0014e-24
9.72768e-25
9.43658e-25
9.14073e-25
8.83997e-25
8.53418e-25
8.22325e-25
7.90708e-25
7.58562e-25
7.25881e-25
6.92664e-25
6.58914e-25
6.24634e-25
5.89833e-25
5.54521e-25
5.18712e-25
4.82424e-25
4.45678e-25
4.08496e-25
3.70906e-25
3.32937e-25
2.94621e-25
2.55992e-25
2.17088e-25
1.77947e-25
1.38611e-25
9.91195e-26
5.95176e-26
1.98488e-26
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