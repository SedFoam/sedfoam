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


/*    This terms compute rho U mean */

time_diff_particle = runtime_particle - T0;

if (time_diff_particle<0)
{
	time_integrate_particle = 0;
}
else
{
	time_integrate_particle = time_diff_particle;
}


timeaverage_particle = time_integrate_particle  + deltaT;


alphaUaMean_a  = (time_integrate_particle* alphaUaMean_a  + alpha* Ua* deltaT) / (timeaverage_particle);

UaMean       = (time_integrate_particle* UaMean + Ua*deltaT) / (timeaverage_particle);

/* alpha Mean */
alpha_aMean_a  = (time_integrate_particle* alpha_aMean_a + alpha* deltaT) / (timeaverage_particle);
/* *****************  Pressure */

pMean_a      = (time_integrate_particle* pMean_a + pa* deltaT ) / (timeaverage_particle);

/* Favre avrage fields */

UaMeanF_a = alphaUaMean_a  / (alpha_aMean_a + 1e-20);

/* Favre fluactuations in velocities and Pressure  */

UaPrimeF_a = Ua - UaMeanF_a;
pPrimeF_a  = pa - pMean_a; 

/* *******************  for UbPrimeMean2F or UaPrimeMean2F     ***********  */

alphaSqrUaMean_a  = (time_integrate_particle* alphaSqrUaMean_a +  alpha* Ua* Ua* deltaT ) / (timeaverage_particle);

UaPrime2MeanF_a = (alphaSqrUaMean_a / (alpha_aMean_a + 1e-20)) - UaMeanF_a* UaMeanF_a;

/* *************** Q and R factor ************ */

volScalarField UaPrimeF_ax = UaPrimeF_a.component(vector::X);
volScalarField UaPrimeF_ay = UaPrimeF_a.component(vector::Y);

RS2_a  = min(UaPrimeF_ax* UaPrimeF_ay* neg(UaPrimeF_ax), zeroo);
RS4_a  = min(UaPrimeF_ax* UaPrimeF_ay* pos(UaPrimeF_ax), zeroo);


RS2Mean_a = (time_integrate_particle* RS2Mean_a + RS2_a* deltaT) / (timeaverage_particle);
RS4Mean_a = (time_integrate_particle* RS4Mean_a + RS4_a* deltaT) / (timeaverage_particle);

/* To write some fundamental variables for TKE budget */


if (runTime.outputTime())
{

	alpha_aMean_a.write();
	UaMeanF_a.write();
	RS2_a.write();
	RS4_a.write();
	UaPrime2MeanF_a.write();
	UaMean.write();
	UaMeanF_a.write();
}

/* ******************* Mean Sij or Tij *********** */


strainTensor_a        = symm(fvc::grad(Ua));
strainTensorM_a       = symm(fvc::grad(UaMeanF_a));
fluctstrainTensor_a   = symm(fvc::grad(UaPrimeF_a));

SijaSija = fluctstrainTensor_a && fluctstrainTensor_a;

SijSij_a   = strainTensorM_a &&  strainTensorM_a;

viscDissMeanI_a = 2* alpha* phasea->nu()* (strainTensor_a && fvc::grad(UaPrimeF_a));
SGSDissMeanI_a  = 2* alpha* turbulencea->nut()* (strainTensor_a && fvc::grad(UaPrimeF_a));


/*           ************** Manohar Formula ************* */

if (favreAverage_particle)
{


	TKEMean_a            =    (time_integrate_particle*TKEMean_a  + deltaT* 0.5* tr(UaPrime2MeanF_a)) / (timeaverage_particle);

	TKEMeanProd_a        =    (time_integrate_particle*TKEMeanProd_a  - deltaT*(UaPrime2MeanF_a && fvc::grad(UaMeanF_a))) / (timeaverage_particle); 

	turbDiffMean_a       =    (time_integrate_particle* turbDiffMean_a - 0.5* alpha* (UaPrimeF_a* magSqr(UaPrimeF_a))* deltaT ) / (timeaverage_particle);
	
	pressDiffMean_a      =    (time_integrate_particle* pressDiffMean_a - (alpha/phasea->rho())* (UaPrimeF_a * p_rbgh)* deltaT ) / (timeaverage_particle); 
	
	pressDiffMeankin_a   =    (time_integrate_particle* pressDiffMeankin_a - (1/phasea->rho())* (UaPrimeF_a * pa)* deltaT ) / (timeaverage_particle); 
	
	pressDiffMeancon_a   =    (time_integrate_particle* pressDiffMeancon_a - (1/phasea->rho())* (UaPrimeF_a * pff)* deltaT ) / (timeaverage_particle); 
	
	viscDiffMean_a       =    (time_integrate_particle* viscDiffMean_a - 2* alpha* (phasea->nu()* (UaPrimeF_a & strainTensor_a ))* deltaT ) / (timeaverage_particle);
	
	viscDissMean_a       =    (time_integrate_particle* viscDissMean_a - (2* alpha* (phasea->nu() + nuFra)* SijaSija)* deltaT) / (timeaverage_particle);
	
//	viscDissMeankin_a    =    (time_integrate_particle* viscDissMeankin_a - (2* alpha* (nua_kin)* SijaSija)* deltaT) / (timeaverage_particle);
	
	viscDissMeancon_a    =    (time_integrate_particle* viscDissMeancon_a - (2* alpha* nuFra* SijaSija)* deltaT) / (timeaverage_particle);
	
	SGSDissMeanJ_a       =    (time_integrate_particle* SGSDissMeanJ_a  - (2* alpha* turbulencea->nut()* SijaSija)* deltaT ) / (timeaverage_particle);

	SGSDiffMean_a        =    (time_integrate_particle* SGSDiffMean_a - 2* alpha* turbulencea->nut()* (UaPrimeF_a & strainTensor_a )* deltaT ) / (timeaverage_particle);

	SGSDissMean_a        =    (time_integrate_particle* SGSDissMean_a  + (2* alpha* turbulencea->nut()* (strainTensor_a && (fvc::grad(UaPrimeF_a)))* deltaT ) )/ (timeaverage_particle);
	
	
	pressDiffMeanJ_a     =    (time_integrate_particle* pressDiffMeanJ_a - (alpha/phasea->rho())* (UaPrimeF_a * pPrimeF_a)* deltaT ) / (timeaverage_particle); 

	SGSDiffMeanJ_a       =    (time_integrate_particle* SGSDiffMeanJ_a - 2* alpha* turbulencea->nut()* (Ua & strainTensor_a )* deltaT ) / (timeaverage_particle);
	
	viscDiffusionMeanJ_a =    (time_integrate_particle* viscDiffusionMeanJ_a +  alpha* phaseb->nu()* (fvc::laplacian(TKEMean_a) )* deltaT ) / (timeaverage_particle);
	
	viscDissMeanJ_a      =    (time_integrate_particle* viscDissMeanJ_a - (2* alpha* phasea->nu()* SijaSija)* deltaT) / (timeaverage_particle);

	SGSDissMeanJ_a       =    (time_integrate_particle* SGSDissMeanJ_a  - (2* alpha* turbulenceb->nut()* SijaSija)* deltaT ) / (timeaverage_particle);

	
}


/* Writing the TKE budget variables */


if (runTime.outputTime())
{

	TKEMean_a.write();
	TKEMeanProd_a.write();
	viscDissMeanI_a.write();
	SGSDissMeanI_a.write();
	turbDiffMean_a.write();
	viscDissMean_a.write();
	viscDissMeanJ_a.write();
//	viscDissMeankin_a.write();
	viscDissMeancon_a.write();
	SGSDissMean_a.write();
	
	turbDiffusionMean_a = fvc::div(turbDiffMean_a);
	pressDiffusionMean_a = fvc::div(pressDiffMean_a);
	pressDiffusionMeankin_a = fvc::div(pressDiffMeankin_a);
	pressDiffusionMeancon_a = fvc::div(pressDiffMeancon_a);
	SGSDiffusionMean_a = fvc::div(SGSDiffMean_a);
	viscDiffusionMean_a = fvc::div(viscDiffMean_a);

        turbDiffusionMean_a.write();
        pressDiffusionMean_a.write();
        SGSDiffusionMean_a.write();
	viscDiffusionMean_a.write();
}


 
/*     *********** Julien formulation  ***************  */

if (favreAverage_particle)
{


	pressDiffMeanJ_a       = (time_integrate_particle* pressDiffMeanJ_a - (alpha/phasea->rho())* (UaPrimeF_a * pPrimeF_a)* deltaT ) / (timeaverage_particle); 

	SGSDiffMeanJ_a         = (time_integrate_particle* SGSDiffMeanJ_a - 2* alpha* turbulencea->nut()* (Ua & strainTensor_a )* deltaT ) / (timeaverage_particle);
	
	viscDiffusionMeanJ_a   = (time_integrate_particle* viscDiffusionMeanJ_a +  alpha* phasea->nu()* (fvc::laplacian(TKEMean_a) )* deltaT ) / (timeaverage_particle);
	
	viscDissMeanJ_a        = (time_integrate_particle* viscDissMeanJ_a - (2* alpha* phasea->nu()* SijaSija)* deltaT) / (timeaverage_particle);

	SGSDissMeanJ_a         = (time_integrate_particle* SGSDissMeanJ_a  - (2* alpha* turbulencea->nut()* SijaSija)* deltaT ) / (timeaverage_particle);
}


/* Writing the TKE budget variables */


if (runTime.outputTime())
{

	viscDissMeanJ_a.write();
	SGSDissMeanJ_a.write();
		
	pressDiffusionMeanJ_a = fvc::div(pressDiffMeanJ_a);
	SGSDiffusionMeanJ_a = fvc::div(SGSDiffMeanJ_a);
	viscDiffusionMeanJ_a = fvc::div(viscDiffMeanJ_a);

        pressDiffusionMeanJ_a.write();
        SGSDiffusionMeanJ_a.write();
	viscDiffusionMeanJ_a.write();
}


runtime_particle +=runTime.deltaTValue();

