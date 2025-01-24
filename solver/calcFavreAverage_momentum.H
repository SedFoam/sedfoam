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

Info<<"Momentum runTime before the addition = \t"<<runtime_momentum <<endl;


time_diff_momentum = runtime_momentum - T0;

//Info<<" Time diff in the favre fluid average file = \t"<<time_diff_momentum<<"\t runtime_momentum  \t"<<runtime_mass<<"\t T0 = \t"<<T0<<endl;

if (time_diff_momentum<0)
{
	time_integrate_momentum = 0;
}
else
{
	time_integrate_momentum = time_diff_momentum;
}


timeaverage_momentum = time_integrate_momentum  + deltaT;

tensor Ei(1.,0,0,0,1.,0,0,0,1.);

Errosion =   (time_integrate_momentum *Errosion + deltaT*alpha* U.component(vector::Y)) / (timeaverage_momentum );
Deposition = (time_integrate_momentum *Deposition +  deltaT* alpha* UfallField) / (timeaverage_momentum );


Taua = (time_integrate_momentum *Taua +  deltaT* fvc::reconstruct(fvc::interpolate(rhoa*(alpha*nuEffa+nuFra)) *mesh.magSf()*fvc::snGrad(Ua)) + rhoa*(alpha*nuEffa+nuFra)*dev(gradUaT)) / (timeaverage_momentum );

Taub = (time_integrate_momentum *Taub +  deltaT* fvc::reconstruct(fvc::interpolate(rhob*(scalar(1.0)-alpha)*nuEffb) *mesh.magSf()*fvc::snGrad(Ub)) + rhob*(scalar(1.0)-alpha)*nuEffb*dev(gradUbT))/(timeaverage_momentum );


volTensorField Rca(
       "Rca",
      ((2.0/3.0)*I)*((nuEffa+nuFra/(alpha+alphaSmall))*tr(gradUaT) + spherSigmaSGSa) - (nuEffa+nuFra/(alpha+alphaSmall))*gradUaT
      );

volTensorField Rca2(
      "Rca2",
       ((2.0/3.0)*I)*((alpha*nuEffa+nuFra)*tr(gradUaT) + alpha*spherSigmaSGSa) - (alpha*nuEffa+nuFra)*gradUaT
       );
        
volVectorField Ra(
       "Ra",
        -nuEffa*fvc::grad(alpha)/(alpha+alphaSmall)
       );
        
surfaceScalarField phiRa(
      -fvc::interpolate(nuEffa)*mesh.magSf()*(fvc::snGrad(alpha))
       /fvc::interpolate(alpha+ alphaSmall)
       );

pMech_b              =    (time_integrate_momentum * pMech_b + deltaT* p)/ (timeaverage_momentum );
pA                   =    (time_integrate_momentum * pA +  deltaT *pa) / (timeaverage_momentum );
pFF                  =    (time_integrate_momentum * pFF + deltaT* pff) / (timeaverage_momentum );
SolidPressureContr_b =    (time_integrate_momentum * SolidPressureContr_b +  deltaT*(-fvc::grad(pa+pff))) / (timeaverage_momentum );
ViscStressContr_b    =   (time_integrate_momentum *ViscStressContr_b + deltaT* rhoa*alpha*fvc::laplacian(nuEffa, Ua)) / (timeaverage_momentum );
FrictContr_b         =  	(time_integrate_momentum * FrictContr_b   +  deltaT* rhoa*fvc::laplacian(nuFra, Ua)) / (timeaverage_momentum );
phiRaContr_b         =    (time_integrate_momentum * phiRaContr_b  +  deltaT* (-rhoa*alpha*fvc::div(phiRa, Ua, "div(phiRa,Ua)"))) / (timeaverage_momentum );
phiRa2Contr_b        =    (time_integrate_momentum * phiRa2Contr_b +  deltaT* rhoa*alpha*fvc::Sp(fvc::div(phiRa), Ua)) / (timeaverage_momentum );
divRcaContr_b        =    (time_integrate_momentum * divRcaContr_b +  deltaT* (-rhoa*alpha*fvc::div(Rca))) / (timeaverage_momentum );
RcaContr_b           =    (time_integrate_momentum * RcaContr_b    +  deltaT* (-rhoa*alpha/(alpha+alphaSmall)*(fvc::grad(alpha) & Rca))) / (timeaverage_momentum );
BulkContr_b          =  	(time_integrate_momentum * BulkContr_b   +  deltaT* fvc::grad(lambdaUa*tr(gradUaT))) / (timeaverage_momentum );
divTauaContr_b       =  ViscStressContr_b + FrictContr_b + phiRaContr_b + phiRa2Contr_b + divRcaContr_b + RcaContr_b + BulkContr_b;
        
gravityContr_b       =  (time_integrate_momentum * gravityContr_b  +    deltaT* alpha/(alpha+alphaSmall)*(rhoa-rhob)*g) / (timeaverage_momentum );
DragContr_b          =  (time_integrate_momentum * DragContr_b     + deltaT* alpha*beta*K*(Ub-Ua)) / (timeaverage_momentum );
ExtContr_b           =  (time_integrate_momentum * ExtContr_b      + deltaT* rhoa*alpha*gradPOSC/(tilt*rhob +(1.0-tilt)*rhoa)) / (timeaverage_momentum );
SuspContr_b          =  (time_integrate_momentum * SuspContr_b     + deltaT* (-SUS*K*beta*turbulenceb->nut()*((SUS_I*iso-SUS_A*aniso) & fvc::grad(alpha)))) / (timeaverage_momentum );

rho_mixd           = (time_integrate_momentum * rho_mixd  + deltaT* (rhoa* alpha* Ua + rhob* beta* Ub)) / (timeaverage_momentum );
nonl_mixd          = (time_integrate_momentum * nonl_mixd + deltaT* (rhoa* alpha* Ua * Ua  + rhob* beta* Ub * Ub)) / (timeaverage_momentum );
bouy_mixd          = (time_integrate_momentum * bouy_mixd + deltaT* (rhoa* alpha  + rhob* beta)* g) / (timeaverage_momentum );

if (runTime.outputTime())
{

Taua.write();
Taub.write();

Errosion.write();
Deposition.write();

Rca.write();
Rca2.write();
Ra.write();


pMech_b.write();
pA.write();
pFF.write();
SolidPressureContr_b.write();
ViscStressContr_b.write();
FrictContr_b.write();
phiRaContr_b.write();
phiRa2Contr_b.write();
divRcaContr_b.write();
RcaContr_b.write();
BulkContr_b.write();
divTauaContr_b.write();
gravityContr_b.write();
DragContr_b.write();
ExtContr_b.write();
SuspContr_b.write();
rho_mixd.write();
nonl_mixd.write();
bouy_mixd.write();
rho.write();

}

runtime_momentum +=runTime.deltaTValue();

