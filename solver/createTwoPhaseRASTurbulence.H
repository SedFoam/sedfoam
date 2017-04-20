dimensionedScalar SUS
(
    twophaseRASProperties.lookup("SUS")
);
dimensionedScalar KE1
(
    twophaseRASProperties.lookup("KE1")
);  
dimensionedScalar Tpsmall
(
    twophaseRASProperties.lookup("Tpsmall")
);  
dimensionedScalar KE3
(
    twophaseRASProperties.lookup("KE3")
);
dimensionedScalar Cmu
(
    twophaseRASProperties.lookup("Cmu")
);
dimensionedScalar B
(
    twophaseRASProperties.lookup("B")
);
volScalarField tmfexp
(
    IOobject
    (
        "tmfexp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -B*rhoa*6.0*Cmu/(max(beta*K,Tpsmall))
);

volScalarField rtemp = (Ub-Ua) & g;
volVectorField Urtemp = rtemp*g/magSqr(g);
volVectorField Urh = (Ub-Ua) - Urtemp;
volScalarField ESD1 = (K/rhob)*(fvc::grad(alpha)) & Urh;

volScalarField ESD3
(
    IOobject
    (
        "ESD3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -2.0*K*alpha/rhob
);

volScalarField ESD4
(
    IOobject
    (
        "ESD4",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (-(SUS/beta)*(rhoa/rhob-1.0)*g & fvc::grad(alpha)) + KE1*ESD1
);

volScalarField ESD5
(
    IOobject
    (
        "ESD5",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    max((-(SUS/beta)*(rhoa/rhob-1.0)*g & fvc::grad(alpha)),scalar(0.0)*ESD1)
    + KE1*ESD1
);

volScalarField ESD
(
    IOobject
    (
        "ESD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    KE3*ESD3
);
