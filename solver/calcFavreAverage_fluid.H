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
    turbDiffusionMean_b    = (time_integrate_fluid* turbDiffusionMean_b    + (fvc::div(turbDiffMean_b)) * (deltaT)) / (timeaverage_fluid);
    turbDiffusionMeanDavid = (time_integrate_fluid* turbDiffusionMeanDavid + (fvc::grad(turbDiffMeanDavid)) * (deltaT) ) / (timeaverage_fluid);
    pressDiffusionMean_b  = (time_integrate_fluid* pressDiffusionMean_b + (fvc::div(pressDiffMean_b))* (deltaT)) / (timeaverage_fluid);
    SGSDiffusionMean_b   = (time_integrate_fluid* SGSDiffusionMean_b + (fvc::div(SGSDiffMean_b))* (deltaT)) / (timeaverage_fluid);
    viscDiffusionMean_b  = (time_integrate_fluid* viscDiffusionMean_b + (fvc::div(viscDiffMean_b))* (deltaT)) / (timeaverage_fluid);
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
    pressDiffusionMeanJ_b  = (time_integrate_fluid* pressDiffusionMeanJ_b + fvc::div(pressDiffMeanJ_b)* deltaT) / (timeaverage_fluid);
    SGSDiffusionMeanJ_b    = (time_integrate_fluid* SGSDiffusionMeanJ_b + fvc::div(SGSDiffMeanJ_b)* deltaT) / (timeaverage_fluid);
    viscDiffusionMeanJ_b   = (time_integrate_fluid* viscDiffusionMeanJ_b + fvc::div(viscDiffMeanJ_b)* deltaT) / (timeaverage_fluid);
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

runtime_fluid += deltaT;
