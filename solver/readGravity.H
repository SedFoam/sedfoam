#include "readGravitationalAcceleration.H"

Info<< "\nStoring g as gEnd\n";
uniformDimensionedVectorField gEnd (g); // gMax already in use in CourantNo.H

// in newer OF versions, this should be expanded to also read g itself
IOdictionary gravitationalProperties
(
    IOobject
    (
        "gravitationalProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

const scalar gRampTime(gravitationalProperties.lookupOrDefault("gRampTime", 
                                                                scalar(0)));

Info << "rampTime for g set to " << gRampTime << "\n";
