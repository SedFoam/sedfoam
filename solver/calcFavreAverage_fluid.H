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

Info<<"Fluid runTime before the addition = \t"<<runtime_fluid<<endl;


time_diff_fluid = runtime_fluid - T0;


if (time_diff_fluid<0)
{
    time_integrate_fluid = 0;
}
else
{
    time_integrate_fluid = time_diff_fluid;
}


deltaT = runTime.deltaTValue();


timeaverage_fluid = time_integrate_fluid  + deltaT;

Info<<"After timeaverage_fluid = \t"<<timeaverage_fluid<<"\t time_integrate_fluid = \t"<<time_integrate_fluid<<endl;
Info<<"After time_diff_fluid = \t"<<time_diff_fluid<<"\t runtime fluid = \t"<<runtime_fluid<<endl;

betaUbMean_F  = (time_integrate_fluid* betaUbMean_F  + beta* Ub* deltaT)  \
                 / (timeaverage_fluid);

UbMean = (time_integrate_fluid* UbMean + Ub*deltaT) / (timeaverage_fluid);

/* alpha Mean */
alpha_aMean_F = (time_integrate_fluid* alpha_aMean_F + alpha* deltaT) \
                / (timeaverage_fluid);

beta_bMean  = (time_integrate_fluid* beta_bMean  + beta* deltaT) \
               / (timeaverage_fluid);
/* *****************  Pressure */

pMean_b = (time_integrate_fluid* pMean_b + p_rbgh* deltaT ) \
              / (timeaverage_fluid);

/* Favre avrage fields */

UbMeanF_b = betaUbMean_F  / (1 - alpha_aMean_F);

/* Favre fluactuations in velocities and Pressure  */

UbPrimeF_b = Ub - UbMeanF_b;
pPrimeF_b  = p_rbgh - pMean_b;

/* *******************  for UbPrimeMean2F or UaPrimeMean2F     ***********  */

betaSqrUbMean_b  = (time_integrate_fluid* betaSqrUbMean_b +  beta* Ub* Ub* deltaT ) \
               / (timeaverage_fluid);

UbPrime2MeanF_b = (betaSqrUbMean_b / (1 - alpha_aMean_F)) \
              -UbMeanF_b* UbMeanF_b;

/* *************** Q and R factor ************ */


volScalarField UbPrimeF_bx = UbPrimeF_b.component(vector::X);
volScalarField UbPrimeF_by = UbPrimeF_b.component(vector::Y);

UbsquarWbMean = (time_integrate_fluid* UbsquarWbMean + (UbPrimeF_bx*UbPrimeF_bx*UbPrimeF_by)* deltaT) / (timeaverage_fluid);
WbcubeMean    = (time_integrate_fluid* WbcubeMean    + (UbPrimeF_by*UbPrimeF_by*UbPrimeF_by)* deltaT) / (timeaverage_fluid);

RS2_b  = min(UbPrimeF_bx* UbPrimeF_by* neg(UbPrimeF_bx), zeroo);
RS4_b  = min(UbPrimeF_bx* UbPrimeF_by* pos(UbPrimeF_bx), zeroo);


RS2Mean_b = (time_integrate_fluid* RS2Mean_b + RS2_b* deltaT) / (timeaverage_fluid);
RS4Mean_b = (time_integrate_fluid* RS4Mean_b + RS4_b* deltaT) / (timeaverage_fluid);

/* To write some fundamental variables for TKE budget */

if (runTime.outputTime())
{
    alpha_aMean_F.write();
    UbMeanF_b.write();
    RS2_b.write();
    RS4_b.write();
    UbPrime2MeanF_b.write();
}

/* ******************* Mean Sij or Tij *********** */

strainTensor_b        = symm(fvc::grad(Ub));
strainTensorM_b = symm(fvc::grad(UbMeanF_b));
fluctstrainTensor_b  = symm(fvc::grad(UbPrimeF_b));

SijbSijb = fluctstrainTensor_b && fluctstrainTensor_b;

SijSij_b   = strainTensorM_b &&  strainTensorM_b;

viscDissMeanI_b = 2* beta* phaseb->nu()* (strainTensor_b && fvc::grad(UbPrimeF_b));
SGSDissMeanI_b  = 2* beta* turbulenceb->nut()* (strainTensor_b && fvc::grad(UbPrimeF_b));

/*           ************** Manohar Formula ************* */

if (favreAverage_fluid)
{
    TKEMeanProdI_b =  - UbPrime2MeanF_b && fvc::grad(UbMeanF_b);
    TKEMean_b         = (time_integrate_fluid* TKEMean_b + deltaT*0.5* tr(UbPrime2MeanF_b)) / (timeaverage_fluid);
    TKEMeanProd_b     = (time_integrate_fluid*TKEMeanProd_b - deltaT* (UbPrime2MeanF_b && fvc::grad(UbMeanF_b))) / (timeaverage_fluid);
    turbDiffMean_b    = (time_integrate_fluid* turbDiffMean_b - 0.5* beta* (UbPrimeF_b* magSqr(UbPrimeF_b))* deltaT ) / (timeaverage_fluid);
    turbDiffMeanDavid = (time_integrate_fluid* turbDiffMeanDavid - 0.5* beta* (UbsquarWbMean + WbcubeMean)* deltaT ) / (timeaverage_fluid);
    pressDiffMean_b   = (time_integrate_fluid* pressDiffMean_b - (beta/phaseb->rho())* (UbPrimeF_b * p_rbgh)* deltaT ) / (timeaverage_fluid);
    SGSDiffMean_b     = (time_integrate_fluid* SGSDiffMean_b - 2* beta* turbulenceb->nut()* (UbPrimeF_b & strainTensor_b )* deltaT ) / (timeaverage_fluid);
    viscDiffMean_b    = (time_integrate_fluid* viscDiffMean_b - 2* beta* phaseb->nu()* (UbPrimeF_b & strainTensor_b )* deltaT ) / (timeaverage_fluid);
    viscDissMean_b    = (time_integrate_fluid* viscDissMean_b - (2* beta* phaseb->nu()* (strainTensor_b && (fvc::grad( UbPrimeF_b))))* deltaT) / (timeaverage_fluid);
    SGSDissMean_b     = (time_integrate_fluid* SGSDissMean_b  + (2* beta* turbulenceb->nut()* (strainTensor_b && (fvc::grad(UbPrimeF_b)))* deltaT ) )/ (timeaverage_fluid);
    dragConstbMean    = (time_integrate_fluid* dragConstbMean - (UbPrimeF_b & (alpha* beta* draga->K(mag(Ub-Ua))* (Ub-Ua))/rhob)* deltaT ) / (timeaverage_fluid);
 
    turbDiffusionMean_b    = (time_integrate_fluid* turbDiffusionMean_b    + (fvc::div(- 0.5* beta* (UbPrimeF_b* magSqr(UbPrimeF_b)))) * (deltaT)) / (timeaverage_fluid);
    turbDiffusionMeanDavid = (time_integrate_fluid* turbDiffusionMeanDavid + (fvc::grad(- 0.5* beta* (UbsquarWbMean + WbcubeMean))) * (deltaT) ) / (timeaverage_fluid);
    pressDiffusionMean_b  = (time_integrate_fluid* pressDiffusionMean_b + (fvc::div(- (beta/phaseb->rho())* (UbPrimeF_b * p_rbgh)))* (deltaT)) / (timeaverage_fluid);
    SGSDiffusionMean_b   = (time_integrate_fluid* SGSDiffusionMean_b + (fvc::div(- 2* beta* turbulenceb->nut()* (UbPrimeF_b & strainTensor_b )))* (deltaT)) / (timeaverage_fluid);
    viscDiffusionMean_b  = (time_integrate_fluid* viscDiffusionMean_b + (fvc::div(- 2* beta* phaseb->nu()* (UbPrimeF_b & strainTensor_b )))* (deltaT)) / (timeaverage_fluid);
   
    //turbDiffusionMean_b    = (time_integrate_fluid* turbDiffusionMean_b    + (fvc::div(turbDiffMean_b)) * (deltaT)) / (timeaverage_fluid);
    //turbDiffusionMeanDavid = (time_integrate_fluid* turbDiffusionMeanDavid + (fvc::grad(turbDiffMeanDavid)) * (deltaT) ) / (timeaverage_fluid);
    //pressDiffusionMean_b  = (time_integrate_fluid* pressDiffusionMean_b + (fvc::div(pressDiffMean_b))* (deltaT)) / (timeaverage_fluid);
    //SGSDiffusionMean_b   = (time_integrate_fluid* SGSDiffusionMean_b + (fvc::div(SGSDiffMean_b))* (deltaT)) / (timeaverage_fluid);
    //viscDiffusionMean_b  = (time_integrate_fluid* viscDiffusionMean_b + (fvc::div(viscDiffMean_b))* (deltaT)) / (timeaverage_fluid);
}

/* Writing the TKE budget variables */

if (runTime.outputTime())

{
    TKEMean_b.write();
    TKEMeanProd_b.write();
    viscDissMeanI_b.write();
    SGSDissMeanI_b.write();
    turbDiffMean_b.write();
    viscDissMean_b.write();
    SGSDissMean_b.write();
    dragConstbMean.write();
    TKEMeanProdI_b.write();
    
    //turbDiffusionMean_b = fvc::div(turbDiffMean_b);
    //turbDiffusionMeanDavid = fvc::grad(turbDiffMeanDavid);
    //pressDiffusionMean_b = fvc::div(pressDiffMean_b);
    //SGSDiffusionMean_b = fvc::div(SGSDiffMean_b);
    //viscDiffusionMean_b = fvc::div(viscDiffMean_b);
    
    turbDiffusionMean_b.write();
    turbDiffusionMeanDavid.write();
    pressDiffusionMean_b.write();
    SGSDiffusionMean_b.write();
    viscDiffusionMean_b.write();

}

/*     *********** Julien formulation  ***************  */

if (favreAverage_fluid)
{

    pressDiffMeanJ_b       = (time_integrate_fluid* pressDiffMeanJ_b - (beta/phaseb->rho())* (UbPrimeF_b * pPrimeF_b)* deltaT ) / (timeaverage_fluid);
    SGSDiffMeanJ_b         = (time_integrate_fluid* SGSDiffMeanJ_b - 2* beta* turbulenceb->nut()* (Ub & strainTensor_b )* deltaT ) / (timeaverage_fluid);
    viscDiffusionMeanJ_b   = (time_integrate_fluid* viscDiffusionMeanJ_b +  beta* phaseb->nu()* (fvc::laplacian(TKEMean_b) )* deltaT ) / (timeaverage_fluid);
    viscDissMeanJ_b        = (time_integrate_fluid* viscDissMeanJ_b - (2* beta* phaseb->nu()* SijbSijb)* deltaT) / (timeaverage_fluid);
    SGSDissMeanJ_b         = (time_integrate_fluid* SGSDissMeanJ_b  - (2* beta* turbulenceb->nut()* SijbSijb)* deltaT ) / (timeaverage_fluid);
 
    pressDiffusionMeanJ_b  = (time_integrate_fluid* pressDiffusionMeanJ_b + fvc::div(- (beta/phaseb->rho())* (UbPrimeF_b * pPrimeF_b))* deltaT) / (timeaverage_fluid);
    SGSDiffusionMeanJ_b    = (time_integrate_fluid* SGSDiffusionMeanJ_b + fvc::div(- 2* beta* turbulenceb->nut()* (Ub & strainTensor_b ))* deltaT) / (timeaverage_fluid);
    //viscDiffusionMeanJ_b   = (time_integrate_fluid* viscDiffusionMeanJ_b + fvc::div(beta* phaseb->nu()* (fvc::laplacian(TKEMean_b) ))* deltaT) / (timeaverage_fluid);
 
    
    //pressDiffusionMeanJ_b  = (time_integrate_fluid* pressDiffusionMeanJ_b + fvc::div(pressDiffMeanJ_b)* deltaT) / (timeaverage_fluid);
    //SGSDiffusionMeanJ_b    = (time_integrate_fluid* SGSDiffusionMeanJ_b + fvc::div(SGSDiffMeanJ_b)* deltaT) / (timeaverage_fluid);
    //viscDiffusionMeanJ_b   = (time_integrate_fluid* viscDiffusionMeanJ_b + fvc::div(viscDiffMeanJ_b)* deltaT) / (timeaverage_fluid);
}

/* Writing the TKE budget variables */

if (runTime.outputTime())

{
    viscDissMeanJ_b.write();
    SGSDissMeanJ_b.write();
    //pressDiffusionMeanJ_b = fvc::div(pressDiffMeanJ_b);
    //SGSDiffusionMeanJ_b = fvc::div(SGSDiffMeanJ_b);
   // viscDiffusionMeanJ_b  = fvc::div(viscDiffMeanJ_b);

    pressDiffusionMeanJ_b.write();
    SGSDiffusionMeanJ_b.write();
    viscDiffusionMeanJ_b.write();

}

/*  Calculation for Mass budget */

if (favreAverage_fluid)
{

alphaF = alpha - alpha_aMean;
Mean_alphaUa     = (time_integrate_fluid* Mean_alphaUa  + alpha* Ua* deltaT) / (timeaverage_fluid);
Mean_UaFalphaF   = (time_integrate_fluid* Mean_UaFalphaF  + alphaF* UaPrimeF_a* deltaT) / (timeaverage_fluid);
alphaFMean       = (time_integrate_fluid* alphaFMean  + alphaF* deltaT) / (timeaverage_fluid);
MeanUaF          = (time_integrate_fluid* MeanUaF  + UaPrimeF_a* deltaT) / (timeaverage_fluid);
}

if (runTime.outputTime())
{
    alphaF.write();
    alphaFMean.write();
    Mean_alphaUa.write();
    Mean_UaFalphaF.write();
    MeanUaF.write();
}


/*  Calculation for particle phase TKE  budget          */


alphaUaMean_a  = (time_integrate_fluid* alphaUaMean_a  + alpha* Ua* deltaT) / (timeaverage_fluid);
UaMean       = (time_integrate_fluid* UaMean + Ua*deltaT) / (timeaverage_fluid);

/* alpha Mean */

alpha_aMean_a  = (time_integrate_fluid* alpha_aMean_a + alpha* deltaT) / (timeaverage_fluid);

/* *****************  Pressure */

pMean_a      = (time_integrate_fluid* pMean_a + pa* deltaT ) / (timeaverage_fluid);

/* Favre avrage fields */

UaMeanF_a = alphaUaMean_a  / (alpha_aMean_a + 1e-20);

/* Favre fluactuations in velocities and Pressure  */

UaPrimeF_a = Ua - UaMeanF_a;
pPrimeF_a  = pa - pMean_a;

/* *******************  for UbPrimeMean2F or UaPrimeMean2F     ***********  */

alphaSqrUaMean_a  = (time_integrate_fluid* alphaSqrUaMean_a +  alpha* Ua* Ua* deltaT ) / (timeaverage_fluid);
UaPrime2MeanF_a = (alphaSqrUaMean_a / (alpha_aMean_a + 1e-20)) - UaMeanF_a* UaMeanF_a;

/* *************** Q and R factor ************ */

volScalarField UaPrimeF_ax = UaPrimeF_a.component(vector::X);
volScalarField UaPrimeF_ay = UaPrimeF_a.component(vector::Y);

RS2_a  = min(UaPrimeF_ax* UaPrimeF_ay* neg(UaPrimeF_ax), zeroo);
RS4_a  = min(UaPrimeF_ax* UaPrimeF_ay* pos(UaPrimeF_ax), zeroo);

RS2Mean_a = (time_integrate_fluid* RS2Mean_a + RS2_a* deltaT) / (timeaverage_fluid);
RS4Mean_a = (time_integrate_fluid* RS4Mean_a + RS4_a* deltaT) / (timeaverage_fluid);

/* To write some fundamental variables for TKE budget */

if (runTime.outputTime())
{
    alpha_aMean_a.write();
    UaMeanF_a.write();
    RS2_a.write();
    RS4_a.write();
    UaPrime2MeanF_a.write();
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

if (favreAverage_fluid)
{
    TKEMeanProdI_a =  - UaPrime2MeanF_a && fvc::grad(UaMeanF_a);
    TKEMean_a            =    (time_integrate_fluid*TKEMean_a  + deltaT* 0.5* tr(UaPrime2MeanF_a)) / (timeaverage_fluid);
    TKEMeanProd_a        =    (time_integrate_fluid*TKEMeanProd_a  - deltaT*(UaPrime2MeanF_a && fvc::grad(UaMeanF_a))) / (timeaverage_fluid);
    
    turbDiffMean_a       =    (time_integrate_fluid* turbDiffMean_a - 0.5* alpha* (UaPrimeF_a* magSqr(UaPrimeF_a))* deltaT ) / (timeaverage_fluid);
    pressDiffMean_a      =    (time_integrate_fluid* pressDiffMean_a - (alpha/phasea->rho())* (UaPrimeF_a * p_rbgh)* deltaT ) / (timeaverage_fluid);
    pressDiffMeankin_a   =    (time_integrate_fluid* pressDiffMeankin_a - (1/phasea->rho())* (UaPrimeF_a * pa)* deltaT ) / (timeaverage_fluid);
    pressDiffMeancon_a   =    (time_integrate_fluid* pressDiffMeancon_a - (1/phasea->rho())* (UaPrimeF_a * pff)* deltaT ) / (timeaverage_fluid);
    pressDiffMeanJ_a     =    (time_integrate_fluid* pressDiffMeanJ_a - (alpha/phasea->rho())* (UaPrimeF_a * pPrimeF_a)* deltaT ) / (timeaverage_fluid);
    SGSDiffMeanJ_a       =    (time_integrate_fluid* SGSDiffMeanJ_a - 2* alpha* turbulencea->nut()* (Ua & strainTensor_a )* deltaT ) / (timeaverage_fluid);     
    viscDiffMean_a       =    (time_integrate_fluid* viscDiffMean_a - 2* alpha* (phasea->nu()* (UaPrimeF_a & strainTensor_a ))* deltaT ) / (timeaverage_fluid);
    SGSDiffMean_a        =    (time_integrate_fluid* SGSDiffMean_a - 2* alpha* turbulencea->nut()* (UaPrimeF_a & strainTensor_a )* deltaT ) / (timeaverage_fluid);    
    
    
    viscDissMean_a       =    (time_integrate_fluid* viscDissMean_a - (2* alpha* (phasea->nu() + nuFra)* SijaSija)* deltaT) / (timeaverage_fluid);
    viscDissMeancon_a    =    (time_integrate_fluid* viscDissMeancon_a - (2* alpha* nuFra* SijaSija)* deltaT) / (timeaverage_fluid);   
    SGSDissMeanJ_a       =    (time_integrate_fluid* SGSDissMeanJ_a  - (2* alpha* turbulencea->nut()* SijaSija)* deltaT ) / (timeaverage_fluid);    
    SGSDissMean_a        =    (time_integrate_fluid* SGSDissMean_a  + (2* alpha* turbulencea->nut()* (strainTensor_a && (fvc::grad(UaPrimeF_a)))* deltaT ) )/ (timeaverage_fluid);
    viscDissMeanJ_a      =    (time_integrate_fluid* viscDissMeanJ_a - (2* alpha* phasea->nu()* SijaSija)* deltaT) / (timeaverage_fluid);
    SGSDissMeanJ_a       =    (time_integrate_fluid* SGSDissMeanJ_a  - (2* alpha* turbulenceb->nut()* SijaSija)* deltaT ) / (timeaverage_fluid);  
 
    
    viscDiffusionMeanJ_a =    (time_integrate_fluid* viscDiffusionMeanJ_a +  alpha* phaseb->nu()* (fvc::laplacian(TKEMean_a) )* deltaT ) / (timeaverage_fluid);      
    turbDiffusionMean_a  =    (time_integrate_fluid* turbDiffusionMean_a + (fvc::div( - 0.5* alpha* (UaPrimeF_a* magSqr(UaPrimeF_a))))* deltaT ) / (timeaverage_fluid);
    pressDiffusionMean_a =    (time_integrate_fluid* pressDiffusionMean_a + (fvc::div( - (alpha/phasea->rho())* (UaPrimeF_a * p_rbgh)))* deltaT ) / (timeaverage_fluid);
    //pressDiffusionMeankin_a =    (time_integrate_fluid* pressDiffusionMeankin_a + (fvc::div(pressDiffMeankin_a))* deltaT ) / (timeaverage_fluid);
    //pressDiffusionMeancon_a =    (time_integrate_fluid* pressDiffusionMeancon_a + (fvc::div(pressDiffMeancon_a))* deltaT ) / (timeaverage_fluid);
    SGSDiffusionMean_a =    (time_integrate_fluid* SGSDiffusionMean_a + (fvc::div(- 2* alpha* turbulencea->nut()* (UaPrimeF_a & strainTensor_a )))* deltaT ) / (timeaverage_fluid);
    viscDiffusionMean_a =    (time_integrate_fluid* viscDiffusionMean_a + (fvc::div(- 2* alpha* (phasea->nu()* (UaPrimeF_a & strainTensor_a ))))* deltaT ) / (timeaverage_fluid);

}

/* Writing the TKE budget variables */

if (runTime.outputTime())
{
    TKEMean_a.write();
    TKEMeanProdI_a.write();
    TKEMeanProd_a.write();
    viscDissMeanI_a.write();
    SGSDissMeanI_a.write();
    turbDiffMean_a.write();
    viscDissMean_a.write();
    viscDissMeanJ_a.write();
    viscDissMeancon_a.write();
    SGSDissMean_a.write();
    
    //turbDiffusionMean_a = fvc::div(turbDiffMean_a);
    //pressDiffusionMean_a = fvc::div(pressDiffMean_a);
    //pressDiffusionMeankin_a = fvc::div(pressDiffMeankin_a);
    //pressDiffusionMeancon_a = fvc::div(pressDiffMeancon_a);
    //SGSDiffusionMean_a = fvc::div(SGSDiffMean_a);
    //viscDiffusionMean_a = fvc::div(viscDiffMean_a);
    
    turbDiffusionMean_a.write();
    pressDiffusionMean_a.write();
    SGSDiffusionMean_a.write();
    viscDiffusionMean_a.write();
}

/*     *********** Julien formulation  ***************  */

if (favreAverage_fluid)
{
    pressDiffMeanJ_a       = (time_integrate_fluid* pressDiffMeanJ_a - (alpha/phasea->rho())* (UaPrimeF_a * pPrimeF_a)* deltaT ) / (timeaverage_fluid);
    SGSDiffMeanJ_a         = (time_integrate_fluid* SGSDiffMeanJ_a - 2* alpha* turbulencea->nut()* (Ua & strainTensor_a )* deltaT ) / (timeaverage_fluid);
    viscDiffusionMeanJ_a   = (time_integrate_fluid* viscDiffusionMeanJ_a +  alpha* phasea->nu()* (fvc::laplacian(TKEMean_a) )* deltaT ) / (timeaverage_fluid);
    viscDissMeanJ_a        = (time_integrate_fluid* viscDissMeanJ_a - (2* alpha* phasea->nu()* SijaSija)* deltaT) / (timeaverage_fluid);
    SGSDissMeanJ_a         = (time_integrate_fluid* SGSDissMeanJ_a  - (2* alpha* turbulencea->nut()* SijaSija)* deltaT ) / (timeaverage_fluid);
    
    pressDiffusionMeanJ_a         = (time_integrate_fluid* pressDiffusionMeanJ_a  +  (fvc::div(- (alpha/phasea->rho())* (UaPrimeF_a * pPrimeF_a)))* deltaT ) / (timeaverage_fluid);
    SGSDiffusionMeanJ_a         = (time_integrate_fluid* SGSDiffusionMeanJ_a  +  (fvc::div(- 2* alpha* turbulencea->nut()* (Ua & strainTensor_a )))* deltaT ) / (timeaverage_fluid);
    viscDiffusionMeanJ_a        = (time_integrate_fluid* viscDiffusionMeanJ_a  +  (fvc::div(viscDiffMeanJ_a))* deltaT ) / (timeaverage_fluid);
}

/* Writing the TKE budget variables */

if (runTime.outputTime())
{
    viscDissMeanJ_a.write();
    SGSDissMeanJ_a.write();
    
    //pressDiffusionMeanJ_a = fvc::div(pressDiffMeanJ_a);
    //SGSDiffusionMeanJ_a = fvc::div(SGSDiffMeanJ_a);
    //viscDiffusionMeanJ_a = fvc::div(viscDiffMeanJ_a);
    
    pressDiffusionMeanJ_a.write();
    SGSDiffusionMeanJ_a.write();
    viscDiffusionMeanJ_a.write();
}

/*   Calculation for the momentum budget */


tensor Ei(1., 0, 0, 0, 1., 0, 0, 0, 1.);

Errosion =   (time_integrate_fluid*Errosion + deltaT*alpha* U.component(vector::Y)) / (timeaverage_fluid );
Deposition = (time_integrate_fluid*Deposition +  deltaT* alpha* UfallField) / (timeaverage_fluid );

Taua = (time_integrate_fluid*Taua +  deltaT* fvc::reconstruct(fvc::interpolate(rhoa*(alpha*nuEffa+nuFra)) *mesh.magSf()*fvc::snGrad(Ua)) + rhoa*(alpha*nuEffa+nuFra)*dev(gradUaT)) / (timeaverage_fluid );

Taub = (time_integrate_fluid*Taub +  deltaT* fvc::reconstruct(fvc::interpolate(rhob*(scalar(1.0)-alpha)*nuEffb) *mesh.magSf()*fvc::snGrad(Ub)) + rhob*(scalar(1.0)-alpha)*nuEffb*dev(gradUbT))/(timeaverage_fluid );

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

pMech_b              =    (time_integrate_fluid* pMech_b + deltaT* p)/ (timeaverage_fluid );
pA                   =    (time_integrate_fluid* pA +  deltaT *pa) / (timeaverage_fluid );
pFF                  =    (time_integrate_fluid* pFF + deltaT* pff) / (timeaverage_fluid );
SolidPressureContr_b =    (time_integrate_fluid* SolidPressureContr_b +  deltaT*(-fvc::grad(pa+pff))) / (timeaverage_fluid );
ViscStressContr_b    =    (time_integrate_fluid*ViscStressContr_b + deltaT* rhoa*alpha*fvc::laplacian(nuEffa, Ua)) / (timeaverage_fluid );
FrictContr_b         =    (time_integrate_fluid* FrictContr_b   +  deltaT* rhoa*fvc::laplacian(nuFra, Ua)) / (timeaverage_fluid );
phiRaContr_b         =    (time_integrate_fluid* phiRaContr_b  +  deltaT* (-rhoa*alpha*fvc::div(phiRa, Ua, "div(phiRa,Ua)"))) / (timeaverage_fluid );
phiRa2Contr_b        =    (time_integrate_fluid* phiRa2Contr_b +  deltaT* rhoa*alpha*fvc::Sp(fvc::div(phiRa), Ua)) / (timeaverage_fluid );
divRcaContr_b        =    (time_integrate_fluid* divRcaContr_b +  deltaT* (-rhoa*alpha*fvc::div(Rca))) / (timeaverage_fluid );
RcaContr_b           =    (time_integrate_fluid* RcaContr_b    +  deltaT* (-rhoa*alpha/(alpha+alphaSmall)*(fvc::grad(alpha) & Rca))) / (timeaverage_fluid );
BulkContr_b          =    (time_integrate_fluid* BulkContr_b   +  deltaT* fvc::grad(lambdaUa*tr(gradUaT))) / (timeaverage_fluid );
divTauaContr_b       =    ViscStressContr_b + FrictContr_b + phiRaContr_b + phiRa2Contr_b + divRcaContr_b + RcaContr_b + BulkContr_b;
gravityContr_b       =    (time_integrate_fluid* gravityContr_b  +    deltaT* alpha/(alpha+alphaSmall)*(rhoa-rhob)*g) / (timeaverage_fluid );
DragContr_b          =    (time_integrate_fluid* DragContr_b     + deltaT* alpha*beta*K*(Ub-Ua)) / (timeaverage_fluid );
ExtContr_b           =    (time_integrate_fluid* ExtContr_b      + deltaT* rhoa*alpha*gradPOSC/(tilt*rhob +(1.0-tilt)*rhoa)) / (timeaverage_fluid );
SuspContr_b          =    (time_integrate_fluid* SuspContr_b     + deltaT* (-SUS*K*beta*turbulenceb->nut()*((SUS_I*iso-SUS_A*aniso) & fvc::grad(alpha)))) / (timeaverage_fluid );

rho_mixd           =      (time_integrate_fluid* rho_mixd  + deltaT* (rhoa* alpha* Ua + rhob* beta* Ub)) / (timeaverage_fluid );
nonl_mixd          =      (time_integrate_fluid* nonl_mixd + deltaT* (rhoa* alpha* Ua * Ua  + rhob* beta* Ub * Ub)) / (timeaverage_fluid );
bouy_mixd          =      (time_integrate_fluid* bouy_mixd + deltaT* (rhoa* alpha  + rhob* beta)* g) / (timeaverage_fluid );

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


/*   Adding time            */
runtime_fluid += deltaT;
