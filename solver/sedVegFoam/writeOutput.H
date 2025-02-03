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

if (favreAveraging)
{
    #include "calcFavreAverage.H"
}

if (favreAverage_fluid)
{
    #include "calcFavreAverage_fluid.H"
    #include "Random.H"
}


if (favreAverage_fluid)
{
    dragVegMean       = (time_integrate_fluid* dragVegMean - (UbPrimeF_b & (beta* alphaVeg* Kveg* (Ub))/rhob)* deltaT ) / (timeaverage_fluid);
}

if (runTime.outputTime())
{

    dragVegMean.write();
}

if (runTime.outputTime())
{

    if (writeTau == 1)
    {
        Taua = fvc::reconstruct(
                    fvc::interpolate(rhoa*(alpha*nuEffa+nuFra)) *
                    mesh.magSf()*fvc::snGrad(Ua)) +
               rhoa*(alpha*nuEffa+nuFra)*dev(gradUaT);
        Taub = fvc::reconstruct(
                    fvc::interpolate(rhob*(scalar(1.0)-alpha-alphaVeg)*nuEffb) *
                    mesh.magSf()*fvc::snGrad(Ub)) +
               rhob*(scalar(1.0)-alpha-alphaVeg)*nuEffb*dev(gradUbT);
        Taua.write();
        Taub.write();
    }
    if (writeMomentumBudget == 1)
    {

        volTensorField Rca(
            "Rca",
            ((2.0/3.0)*I)*((nuEffa+nuFra/(alpha+alphaSmall))*tr(gradUaT)
            + spherSigmaSGSa)
            - (nuEffa+nuFra/(alpha+alphaSmall))*gradUaT
        );

        volTensorField Rca2(
            "Rca2",
            ((2.0/3.0)*I)*((alpha*nuEffa+nuFra)*tr(gradUaT)
            + alpha*spherSigmaSGSa)
            - (alpha*nuEffa+nuFra)*gradUaT
        );
        volVectorField Ra(
            "Ra",
            -nuEffa*fvc::grad(alpha)/(alpha+alphaSmall)
        );

        surfaceScalarField phiRa(
            -fvc::interpolate(nuEffa)*mesh.magSf()*(fvc::snGrad(alpha))
            /fvc::interpolate(alpha+ alphaSmall)
        );

        SolidPressureContr = -fvc::grad(pa+pff);
        ViscStressContr = rhoa*alpha*fvc::laplacian(nuEffa, Ua);
        FrictContr = rhoa*fvc::laplacian(nuFra, Ua);
        phiRaContr = -rhoa*alpha*fvc::div(phiRa, Ua, "div(phiRa,Ua)");
        phiRa2Contr = rhoa*alpha*fvc::Sp(fvc::div(phiRa), Ua);
        divRcaContr = -rhoa*alpha*fvc::div(Rca);
        RcaContr = -rhoa*alpha/(alpha+alphaSmall)*(fvc::grad(alpha) & Rca);
        BulkContr = fvc::grad(lambdaUa*tr(gradUaT));
        divTauaContr = ViscStressContr + FrictContr + phiRaContr +
                       phiRa2Contr + divRcaContr + RcaContr + BulkContr;

        gravityContr=alpha/(alpha+alphaSmall)*(rhoa-rhob)*g;
        DragContr=alpha*beta*K*(Ub-Ua);
        VegDragContr=beta*alphaVeg*Kveg*Ub;
        ExtContr=rhoa*alpha*gradPOSC/(tilt*rhob +(1.0-tilt)*rhoa);
        SuspContr=-SUS*K*beta*turbulenceb->nut()
                  *((SUS_I*iso-SUS_A*aniso) & fvc::grad(alpha));

        divTauaContr.write();
        gravityContr.write();
        DragContr.write();
        VegDragContr.write();
        ExtContr.write();
        SuspContr.write();
    }
    if (granularRheology.on())
    {
        pa.write();
        granularRheology.muI().write();
        if (granularRheology.dilatancy())
        {
            alphaPlastic.write();
            delta.write();
        }
    }
    if (kineticTheory.on())
    {
        pa.write();
        kineticTheory.Theta().write();
    }
    if (SUSlocal)
    {
        SUS.write();
    }
    pff.write();
    runTime.write();
}
