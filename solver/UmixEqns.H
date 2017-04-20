scalar pi = constant::mathematical::pi;

// following is general form of mean flow + sine wave and 2nd order stokes waves
gradPOSC = gradPMEAN+gradPAMP1*Foam::cos(initTheta*(pi)/180.0 + 2.0*(pi)*runTime.value()/oscpT)+gradPAMP2*Foam::sin(initTheta*(pi)/90.0 + 4.0*(pi)*runTime.value()/oscpT);     // Add oscilitary pressure gradient 

//rho=alpha*rhoa+beta*rhob;
rho=alpha*rhoa;

// phase a only rho=rhoa;

volVectorField gradp
    (
        IOobject
        (
            "gradp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        //dimensionedVector("zero", dimensionSet(0, 1, -2, 0, 0), vector::zero)
	dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0), vector::zero)
    );
gradp=gradPMEAN*beta/(beta);

//fvVectorMatrix UaEqn(Ua, Ua.dimensions()*dimVol/dimTime);
fvVectorMatrix UaEqn(Ua, Ua.dimensions()*rho.dimensions()*dimVol/dimTime);

// Effective viscosity of the mixture
nuEffa = nua;
nuEffb = nub;//granularRheology.nuvb();
turbulenceb->nut()   = nuEffb*0.;

// Granular viscosity 
nuFra = granularRheology.mua()/rhoa;


{
    {
        UaEqn =
        (
	 rho*fvm::ddt(Ua)
	 - fvm::laplacian((1-alpha)*(nuEffb+turbulenceb->nut())*rhob, Ua)
	 - fvm::laplacian((nuFra+nua)*rhoa, Ua)
        );
          
        UaEqn.relax();
    }
}

solve(UaEqn == alpha*gradp);
//solve(UaEqn == rho*g);

/* // phase a only 
{
    {
        UaEqn =
        (
         rho*fvm::ddt(Ua)
         -(scalar(1.0)/(alpha+alphaSmall))*fvm::laplacian((nuFra+nua)*rhoa, Ua)
        );

        UaEqn.relax();
    }

}
solve(UaEqn == rho*g);
*/

Ua.correctBoundaryConditions();
