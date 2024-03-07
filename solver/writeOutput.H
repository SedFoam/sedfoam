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
if (runTime.outputTime())
{
    if (writeTau == 1)
    {
        Taua = rhoa*(alpha*nuEffa+nuFra)*dev(fvc::grad(Ua));
        Taub = rhob*(scalar(1.0)-alpha)*nuEffb*dev(fvc::grad(Ub));
        Taua.write();
        Taub.write();
    }
    if (writeMomentumBudget == 1)
    {
        SolidPressureContr=-fvc::grad(pa+pff);
        FluidPressureContr=-fvc::grad(p_rbgh);
        ViscStressContr=rhob*fvc::laplacian(nuEffb, Ub);
        FrictContr=rhoa*fvc::laplacian(nuFra, Ua);
        gravityContr=alpha/(alpha+alphaSmall)*(rhoa-rhob)*g;
        DragContr=alpha*beta*K*(Ub-Ua);
        ExtContr=gradPOSC*rhoa/(tilt*rhob +(1.0-tilt)*rhoa);
        SuspContr=-SUS*K*beta*turbulenceb->nut()
                   *((SUS_I*iso-SUS_A*aniso) & fvc::grad(alpha));;
        
        SolidPressureContr.write();
        FluidPressureContr.write();
        ViscStressContr.write();
        FrictContr.write();
        gravityContr.write();
        DragContr.write();
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
        if (granularRheology.cohesion())
        {
		CohesionDistrb.write();
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
