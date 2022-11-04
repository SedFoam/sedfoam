/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.4.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      turbulenceProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//simulationType  laminar;
//

simulationType  RAS;

RAS{

RASModel        twophasekEpsilon;
//RASModel        twophasekOmega;

turbulence      on;
printCoeffs     on;
    twophasekEpsilonCoeffs
    {
        C1        1.44;
        C2        1.92;
        C3ep      1.2;
        C4ep      1;
        alphak    1.0;
        alphaEps  1.3;
        Cmu              0.09;
        KE2              1.0; //turb modulation
        KE4              1.0; //density stratification g
        nutMax           5e-3;
    }
}
// ************************************************************************* //
