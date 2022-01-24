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

simulationType  LES;

LES{
     LESModel        fluidDynamicLagrangian;

     turbulence      on;

     printCoeffs     on;

     delta           cubeRootVol;

     filter          simple;

     CybMax              0.3;

     nutbMax         1e-4;

     vanDriestCoeffs
     {
         delta           cubeRootVol;
         cubeRootVolCoeffs
         {
             deltaCoeff      1;
         }

         Aplus           26;
         Cdelta          0.158;
     }
}

// ************************************************************************* //

