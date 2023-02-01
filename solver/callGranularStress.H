/*---------------------------------------------------------------------------*\
Copyright (C) 2015 Cyrille Bonamy, Julien Chauchat, Tian-Jian Hsu
                   and contributors

License
    This file is part of SedFOAM.

    SedFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SedFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with SedFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

gradUaT = fvc::grad(Ua)().T();
gradUbT = fvc::grad(Ub)().T();
volScalarField magDtensor
(
    ::sqrt(2.0)*mag(symm(gradUaT))
);

if (granularRheology.dilatancy())
{
    //// Solving the evolution of plastic volume fraction

    fvScalarMatrix phi_pl_Eqn
    (
     fvm::ddt(alphaPlastic)
     + fvm::div(phia, alphaPlastic, "div(phia,alphaPlastic)")
     - fvm::Sp(fvc::div(phia), alphaPlastic)
     ==
     //-fvc::Sp(delta*magDtensor, alphaPlastic) // explicit term
     -fvm::SuSp(delta*magDtensor, alphaPlastic)
    );
    phi_pl_Eqn.relax();
    phi_pl_Eqn.solve();

    alphaPlastic.max(0.40);
    alphaPlastic.min(0.61);
}
else
{
    alphaPlastic=alphaMinFriction*Unity;
}

// Pff corresponds to "contact" pressure which prevent the granular phase from
// exceeding the maximum volume fraction.
pff = pp->pff(alpha, alphaPlastic, alphaMax, Fr_, eta0_, eta1_);

//
//
//  Kinetic THEORY PART
//
//
if (kineticTheory.on())
{
    dimensionedScalar Dsmall2
    (
        "Dsmall2",
        dimensionSet(0, 0, -2, 0, 0, 0, 0),
        1e-8
    );
    Dsmall2 = sqr(Dsmall);

    // Compute Kinetic Theory including granular temperature solution
    kineticTheory.solve
    (
        gradUaT, turbulenceb->k(), turbulenceb->epsilon(),
        turbulenceb->nut(), B, runTime
    );

    // Compute Frictional viscosity
    volScalarField muEff_f
    (
        pff*Foam::sin(kineticTheory.phi())
       /sqrt(pow(magDtensor, 2) + Dsmall2)
    );

    // the actual expression for nuFra depends on the way this term is
    // discretized in UaEqn
    nuFra = muEff_f/rhoa;

    // Compute nuEffb the total fluid phase viscosity
    nuEffb = nub + turbulenceb->nut();

    // Update solid phase viscosities and collisional pressure
    nuEffa = kineticTheory.mua()/((alpha + alphaSmall)*rhoa)
             + turbulencea->nut();
    lambdaUa = kineticTheory.lambda();
    pa = kineticTheory.pa();

    if (debugInfo)
    {
        Info<< "Contact pressure  pff: Min =" << gMin(pff)
            <<", Max =" << gMax(pff)<<endl;
        Info<< "Collisional press. pa: Min =" << gMin(pa)
            <<", Max =" << gMax(pa)<<endl;
    }
}

//
//
//  GRANULAR RHEOLOGY PART
//
//
else if (granularRheology.on())
{
// Solve granular rheology
    granularRheology.solve(magDtensor, pff, alphaSmall, Dsmall);

// Particulate pressure and viscosity
    pa = granularRheology.pa();

    delta=granularRheology.delta();

    if (debugInfo)
    {
        Info<< "Contact pressure  pff: Min =" << gMin(pff)
            <<", Max =" << gMax(pff)<<endl;
        Info<< "Shear ind. press.  pa: Min =" << gMin(pa)
            <<", Max =" << gMax(pa)<<endl;
    }
// the actual expression for nuFra depends on the way this term is discretized
// in UaEqn
    nuFra = granularRheology.mua()/ rhoa;

    nuEffa = sqr(Ct)*turbulenceb->nut() + nua + turbulencea->nut();
    lambdaUa = nua*rhoa*scalar(0.0);

// Compute nuEffb the total fluid phase viscosity
    nuEffb = turbulenceb->nut() + granularRheology.nuvb();
}
//
//
//  If not using kinetic theory or granular rheology Ct model is used
//
//
else
{
    nuEffa = sqr(Ct)*turbulenceb->nut() + nua;
    nuEffb = turbulenceb->nut() + nub;

    // set pa and nuFra to zero if kineticTheory and granularRheology are off
    pa = pa*scalar(0.0);
    nuFra = nua*scalar(0.0);
    nuEffa = nua + turbulencea->nut();
    lambdaUa = nua*rhoa*scalar(0.0);
}
// Add a numerical viscosity to damp the instabilities close to the outlet
if (spongeLayer)
{
    volScalarField XX(mesh.C().component(vector::X));
    nuEffa +=
    (
        pos(XX-xSmin)*nua*1e3
       *Foam::exp(-(xSmax-XX)/max(XX-xSmin, 1e-10*(xSmax-xSmin)))
    );
    nuEffb +=
    (
        pos(XX-xSmin)*nub*1e3
       *Foam::exp(-(xSmax-XX)/max(XX-xSmin, 1e-10*(xSmax-xSmin)))
);
}
pS = pa+pff;

// Limit viscosities for numerical stability
nuFra.min(nuMax);
nuEffa.min(nuMax);
nuEffb.min(nuMax);
