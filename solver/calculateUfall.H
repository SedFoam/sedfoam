volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("eps", dimensionSet(0, 1, -1, 0, 0), 1e-6)
);
volScalarField err
(
    IOobject
    (
        "err",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("err", dimensionSet(0, 1, -1, 0, 0), 0.1)
);

volScalarField Ufalltmp
(
    IOobject
    (
        "Ufalltmp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("small_", dimensionSet(0, 1, -1, 0, 0), 1.e-6)
);


volScalarField UfallField
(
    IOobject
    (
        "UfallField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimensionSet(0, 1, -1, 0, 0), 1.0)
);

scalar iter=0;
while ((err > eps)||(iter < 10))
{
    iter = iter + 1;
    UfallField = mag(g)*(rhoa-rhob)/draga->K(Ufalltmp);
    err = max(mag(UfallField - Ufalltmp)); 
    // under-relaxation for convergence stability
    Ufalltmp = 0.5*Ufalltmp+0.5*UfallField;
    //Info<<iter<<endl;
}

if (iter == 10)
{
    Info<<"Maximum number of iterations reached for settling velocity calculation"<<endl;
}
Info<<"Settling velocity in quiescent water: Ufall="<<max(UfallField).value()<<" with residual:"<<err<<endl;
