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

Info<<"Compute Favre-averaged quantities"<<endl;

if (runTime.value() < T0)
{
    time_average = 0;
}
else
{
    deltaT = runTime.deltaTValue();
    time_average = runTime.value() - T0 + deltaT;
    
    betaUbMean  = ((time_average - deltaT) * betaUbMean  + deltaT * beta* Ub ) \
                  / time_average;

    /* alpha Mean */
    alpha_aMean = ((time_average - deltaT) * alpha_aMean + deltaT * alpha ) \
                  / time_average;

    /* *****************  Pressure */
    pMean = ((time_average - deltaT) * pMean + deltaT * p_rbgh ) \
            / time_average;

    /* Favre avrage fields */
    UbMeanF = betaUbMean  / (1 - alpha_aMean);

    /* Favre fluactuations in velocities and Pressure  */
    UbPrimeF = Ub - UbMeanF;
    pPrime  = p_rbgh - pMean;

    /* ******** for UbPrimeMean2F or UaPrimeMean2F ********  */
    betaSqrUbMean  = ((time_average - deltaT) * betaSqrUbMean \
                      +  deltaT * beta * Ub * Ub ) \
                     / time_average;

    UbPrime2MeanF = (betaSqrUbMean / (1 - alpha_aMean)) \
                    - UbMeanF * UbMeanF;

    /* ***************  R factor ************ */
    volScalarField UbPrimeF_x = UbPrimeF.component(vector::X);
    volScalarField UbPrimeF_y = UbPrimeF.component(vector::Y);
    
    UbSqrWbMean = ((time_average - deltaT) * UbSqrWbMean \
                   + deltaT * (UbPrimeF_x * UbPrimeF_x * UbPrimeF_y) ) \
                   / time_average;

    WbcubeMean    = ((time_average - deltaT) * WbcubeMean \
                     + deltaT * (UbPrimeF_y * UbPrimeF_y * UbPrimeF_y) )
                    / time_average;

    RS2_b  = min(UbPrimeF_x* UbPrimeF_y* neg(UbPrimeF_x), zeroo);
    RS4_b  = min(UbPrimeF_x* UbPrimeF_y* pos(UbPrimeF_x), zeroo);

    RS2Mean_b = ((time_average - deltaT) * RS2Mean_b \
                 + deltaT * RS2_b ) \
                / time_average;

    RS4Mean_b = ((time_average - deltaT) * RS4Mean_b \
                 + deltaT * RS4_b ) \
                / time_average;

    /* To write some fundamental variables for TKE budget */
    if (runTime.outputTime())
    {
        alpha_aMean.write();
        betaUbMean.write();
        betaSqrUbMean.write();
        UbMeanF.write();
        RS2_b.write();
        RS4_b.write();
        UbPrime2MeanF.write();
    }

    if (TKEbudgetLES_fluid)
    {
        Info<<"Compute Turbulent Kinetic Energy Budget (for LES or DNS)"<<endl;
        /* ******************* Mean Sij or Tij *********** */
        strainTensor_b = symm(fvc::grad(Ub));
        strainTensorM_b = symm(fvc::grad(UbMeanF));
        fluctstrainTensor_b = symm(fvc::grad(UbPrimeF));
        
        SijbSijb = fluctstrainTensor_b && fluctstrainTensor_b;
        SijSij_b = strainTensorM_b &&  strainTensorM_b;
        
        /*           ************** Manohar Formula ************* */
        
        TKEMean_b = ((time_average - deltaT) * TKEMean_b
                    + deltaT*0.5* tr(UbPrime2MeanF))
                    / time_average;

        TKEMeanProd_b = ((time_average - deltaT) * TKEMeanProd_b
                     - deltaT* (UbPrime2MeanF && fvc::grad(UbMeanF)))
                      / time_average;

        dragConstbMean    = ((time_average - deltaT) * dragConstbMean
                        - (UbPrimeF & (alpha* beta* draga->K(mag(Ub-Ua))*
                          (Ub-Ua))/rhob)* deltaT ) / time_average;

        turbDiffusionMean_b = ((time_average - deltaT) * turbDiffusionMean_b
                         + (fvc::div(- 0.5* beta* (UbPrimeF*
                           magSqr(UbPrimeF))))
                         * (deltaT)) / time_average;

        turbDiffusionMeanWallNormal_b = ((time_average - deltaT) *
                                      turbDiffusionMeanWallNormal_b
                         + (fvc::grad(- 0.5* beta*(UbSqrWbMean + WbcubeMean)))
                         * (deltaT)) / time_average;

        viscDiffusionMean_b = ((time_average - deltaT) * viscDiffusionMean_b
                             +  beta* phaseb->nu()* (fvc::laplacian(TKEMean_b))
                             * deltaT ) / time_average;

        viscDissMean_b      = ((time_average - deltaT) * viscDissMean_b
                             - (2* beta* phaseb->nu()* SijbSijb)* deltaT)
                             / time_average;

        SGSDissMean_b       = ((time_average - deltaT) * SGSDissMean_b
                             - (2* beta* turbulenceb->nut()* SijbSijb)* deltaT)
                             / time_average;

        pressDiffusionMean_b = ((time_average - deltaT) * pressDiffusionMean_b
                             + fvc::div(- (beta/phaseb->rho())* (UbPrimeF
                             * pPrime))* deltaT) / time_average;

        SGSDiffusionMean_b   = ((time_average - deltaT) * SGSDiffusionMean_b
                             + fvc::div(- 2* beta* turbulenceb->nut()
                             * (Ub & strainTensor_b ))* deltaT)
                             / time_average;

        /* Writing the TKE budget variables */
        if (runTime.outputTime())
        {
            TKEMean_b.write();
            TKEMeanProd_b.write();
            dragConstbMean.write();
            turbDiffusionMean_b.write();
            turbDiffusionMeanWallNormal_b.write();
             viscDiffusionMean_b.write();
            viscDissMean_b.write();
            SGSDissMean_b.write();
            pressDiffusionMean_b.write();
            SGSDiffusionMean_b.write();
            RS2Mean_b.write();
            RS4Mean_b.write();
        }
    }
    /*  Calculation for Mass budget */
    if (favreAverage_solid)
    {
        /*  Calculation for particle phase TKE  budget          */
        alphaUaMean  = ((time_average - deltaT) * alphaUaMean
                        + deltaT * alpha* Ua )
                       / time_average;

        /* *****************  Pressure */
        paMean       = ((time_average - deltaT) * pMean + deltaT * (pa+pff) )
                       / time_average;

        /* Favre avrage fields */
        UaMeanF = alphaUaMean  / (alpha_aMean + alphaSmall);

        /* Favre fluactuations in velocities and Pressure  */
        UaPrimeF = Ua - UaMeanF;
        paPrime  = pa+pff - paMean;

        /* ***  for UbPrimeMean2F or UaPrimeMean2F     ********  */
        alphaSqrUaMean  = ((time_average - deltaT) * alphaSqrUaMean
                           + deltaT *alpha * Ua * Ua )
                          / time_average;

        UaPrime2MeanF = (alphaSqrUaMean / (alpha_aMean + alphaSmall))
                         - UaMeanF * UaMeanF;

        /* *************** Q and R factor ************ */
        volScalarField UaPrimeF_x = UaPrimeF.component(vector::X);
        volScalarField UaPrimeF_y = UaPrimeF.component(vector::Y);

        RS2_a  = min(UaPrimeF_x * UaPrimeF_y* neg(UaPrimeF_x), zeroo);
        RS4_a  = min(UaPrimeF_x * UaPrimeF_y* pos(UaPrimeF_x), zeroo);

        RS2Mean_a = ((time_average - deltaT) * RS2Mean_a + RS2_a* deltaT)
                    /time_average;
        RS4Mean_a = ((time_average - deltaT) * RS4Mean_a + RS4_a* deltaT)
                    / time_average;

        /* To write some fundamental variables for TKE budget */
        if (runTime.outputTime())
        {
            alphaUaMean.write();
            alphaSqrUaMean.write();
            UaMeanF.write();
            RS2_a.write();
            RS4_a.write();
            UaPrime2MeanF.write();
        }
    }
    if (TKEbudgetLES_solid)
    {
        /* ******************* Mean Sij or Tij *********** */
        strainTensor_a        = symm(fvc::grad(Ua));
        strainTensorM_a       = symm(fvc::grad(UaMeanF));
        fluctstrainTensor_a   = symm(fvc::grad(UaPrimeF));

        SijaSija = fluctstrainTensor_a && fluctstrainTensor_a;
        SijSij_a   = strainTensorM_a &&  strainTensorM_a;

        TKEMean_a      =  ((time_average - deltaT) * TKEMean_a
                              + deltaT* 0.5* tr(UaPrime2MeanF))
                          / time_average;
        TKEMeanProd_a  =  ((time_average - deltaT) * TKEMeanProd_a
                              - deltaT*(UaPrime2MeanF && fvc::grad(UaMeanF)))
                          / time_average;

        turbDiffusionMean_a = ((time_average - deltaT) * turbDiffusionMean_a
                       + (fvc::div(- 0.5*alpha*(UaPrimeF * magSqr(UaPrimeF))))*
                       deltaT ) / time_average;

        viscDiffusionMean_a = ((time_average - deltaT) * viscDiffusionMean_a
                    +  alpha* phaseb->nu()* (fvc::laplacian(TKEMean_a) )*
                       deltaT ) / time_average;

        viscDissMean_a   = ((time_average - deltaT) * viscDissMean_a
                     - (2* alpha* (phasea->nu() + nuFra)* SijaSija)* deltaT)
                     / time_average;

        SGSDissMean_a = ((time_average - deltaT) * SGSDissMean_a
                     - (2* alpha* turbulencea->nut()* SijaSija)* deltaT )
                     / time_average;

        pressDiffusionMean_a  = ((time_average - deltaT) * pressDiffusionMean_a
                             + (fvc::div(- (alpha/phasea->rho())*
                               (UaPrimeF * paPrime)))* deltaT )
                             / time_average;

        SGSDiffusionMean_a    = ((time_average - deltaT) * SGSDiffusionMean_a
                             + (fvc::div(- 2* alpha* turbulencea->nut()*
                               (Ua & strainTensor_a )))* deltaT )
                             / time_average;

        /* Writing the TKE budget variables */
        if (runTime.outputTime())
        {
            TKEMean_a.write();
            TKEMeanProd_a.write();
            turbDiffusionMean_a.write();
            viscDiffusionMean_a.write();
            viscDissMean_a.write();
            SGSDissMean_a.write();
            pressDiffusionMean_a.write();
            SGSDiffusionMean_a.write();
            RS2Mean_a.write();
            RS4Mean_a.write();
        }
    }

/*   Calculation for the momentum budget */
/* Manohar's version of momentume budget - commented for the moment
// unit tensor
tensor Ei(1., 0, 0, 0, 1., 0, 0, 0, 1.);

Taua = ((time_average - deltaT)*Taua +  deltaT*
            fvc::reconstruct(fvc::interpolate(rhoa*(alpha*nuEffa+nuFra))*
            mesh.magSf()*fvc::snGrad(Ua)) + rhoa*(alpha*nuEffa+nuFra)*
            dev(gradUaT)) / (timeaverage_fluid );

Taub = ((time_average - deltaT)*Taub +  deltaT*
            fvc::reconstruct(fvc::interpolate(rhob*(scalar(1.0)-alpha)*nuEffb)*
            mesh.magSf()*fvc::snGrad(Ub)) + rhob*(scalar(1.0)-alpha)*
            nuEffb*dev(gradUbT))/(timeaverage_fluid );

volTensorField Rca(
    "Rca",
    ((2.0/3.0)*I)*((nuEffa+nuFra/(alpha+alphaSmall))*tr(gradUaT)
     + spherSigmaSGSa) - (nuEffa+nuFra/(alpha+alphaSmall))*gradUaT
    );

volTensorField Rca2(
    "Rca2",
    ((2.0/3.0)*I)*((alpha*nuEffa+nuFra)*tr(gradUaT) + alpha*spherSigmaSGSa)
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

pMech_b = ((time_average - deltaT)* pMech_b + deltaT* p)/ (timeaverage_fluid );
pA      = ((time_average - deltaT)* pA +  deltaT *pa) / (timeaverage_fluid );
pFF     = ((time_average - deltaT)* pFF + deltaT* pff) / (timeaverage_fluid );
SolidPressureContr_b =  ((time_average - deltaT)* SolidPressureContr_b
                        +  deltaT*(-fvc::grad(pa+pff))) / (timeaverage_fluid );
ViscStressContr_b = ((time_average - deltaT)*ViscStressContr_b
                    + deltaT* rhoa*alpha*fvc::laplacian(nuEffa, Ua))
                    / (timeaverage_fluid );
FrictContr_b = ((time_average - deltaT)* FrictContr_b   +  deltaT* rhoa*
                fvc::laplacian(nuFra, Ua)) / (timeaverage_fluid );
phiRaContr_b = ((time_average - deltaT)* phiRaContr_b  +  deltaT*
               (-rhoa*alpha*fvc::div(phiRa, Ua, "div(phiRa,Ua)")))
               / (timeaverage_fluid );
phiRa2Contr_b = ((time_average - deltaT)* phiRa2Contr_b +  deltaT* rhoa*
                alpha*fvc::Sp(fvc::div(phiRa), Ua)) / (timeaverage_fluid );
divRcaContr_b = ((time_average - deltaT)* divRcaContr_b +  deltaT*
                (-rhoa*alpha*fvc::div(Rca))) / (timeaverage_fluid );
RcaContr_b = ((time_average - deltaT)* RcaContr_b    +  deltaT*
                (-rhoa*alpha/(alpha+alphaSmall)*(fvc::grad(alpha) & Rca)))
                / (timeaverage_fluid );
BulkContr_b = ((time_average - deltaT)* BulkContr_b   +  deltaT*
                  fvc::grad(lambdaUa*tr(gradUaT))) / (timeaverage_fluid );
divTauaContr_b = ViscStressContr_b + FrictContr_b + phiRaContr_b
                 + phiRa2Contr_b + divRcaContr_b + RcaContr_b + BulkContr_b;
gravityContr_b = ((time_average - deltaT)* gravityContr_b  +    deltaT*
                  alpha/(alpha+alphaSmall)*(rhoa-rhob)*g)
                  / (timeaverage_fluid );
DragContr_b = ((time_average - deltaT)* DragContr_b + deltaT* alpha
               *beta*K*(Ub-Ua)) / (timeaverage_fluid );
ExtContr_b = ((time_average - deltaT)* ExtContr_b + deltaT* rhoa*alpha*
                gradPOSC/(tilt*rhob +(1.0-tilt)*rhoa)) / (timeaverage_fluid );
SuspContr_b = ((time_average - deltaT)* SuspContr_b + deltaT* (-SUS*K*beta*
               turbulenceb->nut()*((SUS_I*iso-SUS_A*aniso) & fvc::grad(alpha))))
               / (timeaverage_fluid );

rho_mixd = ((time_average - deltaT)* rho_mixd  + deltaT* (rhoa* alpha* Ua
               + rhob* beta* Ub)) / (timeaverage_fluid );
nonl_mixd = ((time_average - deltaT)* nonl_mixd + deltaT* (rhoa* alpha* Ua * Ua
               + rhob* beta* Ub * Ub)) / (timeaverage_fluid );
buoy_mixd = ((time_average - deltaT)* buoy_mixd + deltaT* (rhoa* alpha
               + rhob* beta)* g) / (timeaverage_fluid );

if (runTime.outputTime())
{
    Taua.write();
    Taub.write();

    Erosion.write();
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
    buoy_mixd.write();
    rho.write();
}
     end Manohar's version for momentum budget */
}
