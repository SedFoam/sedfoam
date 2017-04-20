volTensorField gradUaT = fvc::grad(Ua)().T();
volTensorField gradUbT = fvc::grad(Ub)().T();
volSymmTensorField devS = dev(symm(gradUaT));
// is galpha really used?
surfaceScalarField galpha =
(
    fvc::interpolate(1.0/(alpha + alphaSmall))*fvc::snGrad(alpha)*mesh.magSf()
);

if (kineticTheory.on())
{
    // Compute Kinetic Theory including granular temperature solution
    kineticTheory.solve
    (
        galpha, gradUaT, turbulenceb->k(), turbulenceb->epsilon(), turbulenceb->nut(),
        B, runTime
    );

    // Update solid phase viscosities and collisional pressure 
    nuEffa = kineticTheory.mua()/((alpha + alphaSmall)*rhoa);
    lambdaUa = kineticTheory.lambda();
    pa = kineticTheory.pa();

    /*    
    // Compute the gradient of collisional pressure (not used in the current version)
    surfaceScalarField gradPsf = fvc::snGrad(pa)*mesh.magSf();
    forAll(p.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            gradPsf.boundaryField()[patchi] = 0.0;
        }
    } 
    gradPs = (fvc::reconstruct(gradPsf))/(rhoa*(alpha+alphaSmall));
    */
}
else               // If not using kinetic theory is using Ct model
{
    nuEffa = sqr(Ct)*turbulenceb->nut() + nua ;
}

if (kineticTheory.on() && debugInfo)
{
    Info<< "Max gradPs =" << max(mag(gradPs)).value()
        << ",Min gradPs =" << min(mag(gradPs)).value() << endl;
    Info<< "Max gradPf =" << max(mag(gradPf)).value()
        << ",Min gradPf =" << min(mag(gradPf)).value() << endl;
}

