//
// Pff corresponds to "contact" pressure which prevent the granular phase from
// exceeding the maximum volume fraction.
// The values of alphaMax, alphaMinFriction, Fr_, eta0_ and eta1_ are read from the ppProperties file
pff =
(
    Fr_*pow(max(alpha - alphaMinFriction, scalar(0)), eta0_)
    /pow(max(alphaMax - alpha, scalar(1.0e-20)), eta1_)
);
//
// Compute the gradient of pff that enter the momentum balance equation
//
surfaceScalarField gradPff = fvc::snGrad(pff)*mesh.magSf();

forAll(p_rbgh.boundaryField(), patchi)
{
    if (isA<zeroGradientFvPatchScalarField>(p_rbgh.boundaryField()[patchi]))
    {
        gradPff.boundaryFieldRef()[patchi] = 0.0;
    }
}

gradPf = (fvc::reconstruct(gradPff))/(rhoa*(alpha + alphaSmall));
//
//
//  Kinetic THEORY PART
//
//
if (kineticTheory.on())
{
    volTensorField dU = fvc::grad(Ua);
    volSymmTensorField dUU = symm(dU);
    volSymmTensorField devS = dUU - (scalar(1.0)/scalar(3.0))*tr(dUU)*I;
    dimensionedScalar I2Dsmall
    (
        "small",
        dimensionSet(0 , 0 ,-2 ,0 , 0, 0, 0),
        scalar(1.0e-40)
    );

    //scalar pi = constant::mathematical::pi;
    // T && T = T_ij*T_ij
    volScalarField muEff_f =
    (
        pff*Foam::sin(kineticTheory.phi())
       /sqrt(scalar(2.0)*(devS && devS) + I2Dsmall)
    );
// the actual expression for nuFra depends on the way this term is discretized
// in UaEqn
    nuFra = muEff_f/rhoa;
//   nuFra = muEff_f/((alpha+alphaSmall)*rhoa);

// Compute nuEffb the total fluid phase viscosity
    nuEffb = nub + turbulenceb->nut();
}

//
//
//  GRANULAR RHEOLOGY PART
//
//
else if (granularRheology.on())
{
    
// solving granular rheology
    granularRheology.solve(gradUaT,pff,alphaSmall,runTime.deltaT());

// Particulate pressure and viscosity
    pa = granularRheology.pa();
    if (debugInfo)
    {
        Info<< "Pff pressure: Max =" << gMax(pff)
            <<",Min =" << gMin(pff)<<endl;
        Info<< "Shear induced pressure : Max =" << gMax(pa)
            <<",Min =" << gMin(pa)<<endl;
    }
/*
// Compute the shear induced pressure gradient contribution (not used inthis version)
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
// the actual expression for nuFra depends on the way this term is discretized
// in UaEqn
    nuFra = granularRheology.mua()/ rhoa;
//     nuFra = granularRheology.mua()/((alpha+alphaSmall)*rhoa);

    nuEffa = sqr(Ct)*turbulenceb->nut() + nua;

// Compute nuEffb the total fluid phase viscosity
    nuEffb = turbulenceb->nut() + granularRheology.nuvb();
}
//
//
//  If not using kinetic theory is using Ct model
//
//
else               
{
    nuEffa = sqr(Ct)*turbulenceb->nut() + nua;
    nuEffb = turbulenceb->nut() + nub;
    
    // set pa and nuFra to zero if kineticTheory and granularRheology are off
    pa = pa*scalar(0.0);
    nuFra = nua*scalar(0.0);
    nuEffa = nua;
}
// Add a numerical viscosity to damp the instabilities close to the outlet
if (spongeLayer)
{
    volScalarField XX = mesh.C().component(vector::X);
    nuEffa +=
    (
        pos(XX-xSmin)*nua*1e3
       *Foam::exp(-(xSmax-XX)/max(XX-xSmin,1e-10*(xSmax-xSmin)))
    );
    nuEffb +=
    (
        pos(XX-xSmin)*nub*1e3
       *Foam::exp(-(xSmax-XX)/max(XX-xSmin,1e-10*(xSmax-xSmin)))
    );
}

// Limit viscosities for numerical stability
nuFra.min(nuMax);
nuEffa.min(nuMax);
nuEffb.min(nuMax);
