IOdictionary ppProperties
(
    IOobject
    (
        "ppProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

scalar alphaMax(readScalar(ppProperties.lookup("alphaMax")));
dimensionedScalar alphaMinFriction(ppProperties.lookup("alphaMinFriction"));
dimensionedScalar Fr_(ppProperties.lookup("Fr"));
dimensionedScalar eta0_
(
    ppProperties.lookupOrDefault
    (
       "eta0_",
       dimensionedScalar("eta0_",dimensionSet(0,0,0,0,0,0,0), 3 )
    )
);
dimensionedScalar eta1_
(
    ppProperties.lookupOrDefault
    (
       "eta1_",
       dimensionedScalar("eta1_",dimensionSet(0,0,0,0,0,0,0), 5 )
    )
);
Switch packingLimiter(ppProperties.lookup("packingLimiter"));
