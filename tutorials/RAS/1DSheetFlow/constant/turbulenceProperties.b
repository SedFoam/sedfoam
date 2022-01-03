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

simulationType  RAS;

RAS{
//RASModel        twophaseMixingLength;
    //RASModel        twophasekEpsilon;
RASModel        twophasekOmega;

    turbulence      on;
    printCoeffs     on;
    twophasekOmegaCoeffs
    {


        alphaOmega       0.52;
        betaOmega        0.072; // original k-omega: 0.072 // k-omega (2006): 0.0708
        C3om             0.35;
        C4om             1.0;
        alphaKomega      0.5; // original k-omega: 0.5 // k-omega (2006): 0.6
        alphaOmegaOmega  0.5;
        Clim             0.0; // original k-omega: 0.0 // k-omega (2006): 0.875
        sigmad           0.0; // original k-omega: 0.0 // k-omega (2006): 0.125
        Cmu              0.09;
        KE2              1.0; //turb modulation
        KE4              1.0; //density stratification g  
        nutMax           5e-3;
        popeCorrection   false;
    }
    twophaseMixingLengthCoeffs
    {
        expoLM       1.0;
        alphaMaxLM   0.55;
        kappaLM      0.225;
        Cmu              0.09;
        nutMax           5e-3;
    }
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
