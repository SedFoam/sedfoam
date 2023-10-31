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

betaUbMean = (beta*Ub + avCount*betaUbMean)/(avCount+1);
alphaUaMean = (alpha*Ua + avCount*alphaUaMean)/(avCount+1);
alpha_aMean = (alpha + avCount*alpha_aMean)/(avCount+1);

pMeanF = (p_rbgh + avCount*pMeanF)/(avCount+1);

UbMeanF = betaUbMean/(1-alpha_aMean);
UaMeanF = alphaUaMean/(alpha_aMean+1e-20);

UbPrimeF = Ub - UbMeanF;
UaPrimeF = Ua - UaMeanF;

pPrimeF = p_rbgh - pMeanF;

volScalarField UbPrimeFx = UbPrimeF.component(vector::X);
volScalarField UbPrimeFy = UbPrimeF.component(vector::Y);

dimensionedScalar zeroo
(
"zeroo",
dimensionSet(0, 2, -2, 0, 0, 0, 0),
scalar(0.0)
);

RS2 = min(UbPrimeFx*UbPrimeFy*neg(UbPrimeFx), zeroo);
RS4 = min(UbPrimeFx*UbPrimeFy*pos(UbPrimeFx), zeroo);

RS2Mean = (RS2 + avCount*RS2Mean)/(avCount+1);
RS4Mean = (RS4 + avCount*RS4Mean)/(avCount+1);

betaSqrUbMean = (beta*Ub*Ub + avCount*betaSqrUbMean)/(avCount+1);
alphaSqrUaMean = (alpha*Ua*Ua + avCount*alphaSqrUaMean)/(avCount+1);

UbPrime2MeanF = (betaSqrUbMean/(1-alpha_aMean))-UbMeanF*UbMeanF;
UaPrime2MeanF = alphaSqrUaMean/(alpha_aMean+1e-20)-UaMeanF*UaMeanF;

if (runTime.outputTime())
{
    alpha_aMean.write();
    RS2Mean.write();
    RS4Mean.write();
    UbMeanF.write();
    UaMeanF.write();
    UbPrime2MeanF.write();
    UaPrime2MeanF.write();
}

if (TKEBudgetLES)
{
    Info<<"Compute TKE budget for LES"<<endl;
    TKEMean = 0.5 * tr(UbPrime2MeanF);
//    Info<<"Compute TKE mean Prod"<<endl;

    TKEMeanProd = - UbPrime2MeanF && fvc::grad(UbMeanF);
//    Info<<"Strain Rate Tensor of Resolved Fluctuations"<<endl;
    fluctstrainTensor = twoSymm(fvc::grad(UbPrimeF))
                      -((2.0/3.0)*I)*tr(fvc::grad(UbPrimeF));
//    Info<<"Strain Rate Tensor of the resolved flow"<<endl;
    strainTensor = twoSymm(fvc::grad(Ub))
                      -((2.0/3.0)*I)*tr(fvc::grad(Ub));

// For the next three terms only the fluxes are time-averaged, 
//   the divergence is taken after averaging
//    Info<<"Turbulent Transport Term"<<endl;
    turbDiffMean = (avCount*turbDiffMean 
                         -0.5*beta*UbPrimeF*magSqr(UbPrimeF)
                        )/(avCount+1);
//    Info<<"Pressure Trasnport Term"<<endl;
    pressDiffMean = (avCount*pressDiffMean 
                          - beta/phaseb->rho()*(UbPrimeF * pPrimeF)
                         )/(avCount+1);
//    Info<<"Subgrid Diffusion Term"<<endl;
    SGSDiffMean = (avCount*SGSDiffMean 
                   - beta * turbulenceb->nut() * (Ub & strainTensor)
                  )/(avCount+1);
//    Info<<"Viscous Diffusion"<<endl; 
    viscDiffusionMean = beta * phaseb->nu()*fvc::laplacian(TKEMean);
//    Info<<"Viscous Dissipation of Resolved Fluctuations"<<endl; 
    viscDissMean = (avCount*viscDissMean 
                    - 0.5 * beta * phaseb->nu()
                          * (fluctstrainTensor && fluctstrainTensor)
                   )/(avCount+1);
//    Info<<"Subgrid Dissipation of Resolved Fluctuations"<<endl; 
    SGSDissMean = (avCount*SGSDissMean 
                   - 0.5*beta * turbulenceb->nut()
                          *(fluctstrainTensor && fluctstrainTensor)
                  )/(avCount+1);
//    Info<<"Turbulent Drag Dissipation/Transfer term"<<endl; 
    dragContbMean = (avCount*dragContbMean
                     - (UbPrimeF & 
                         (alpha*beta*draga->K(mag(Ub-Ua))*(Ub-Ua))/rhob)
                    )/(avCount+1);
//    dragContbMean = ((-UbPrimeF & (alpha*beta*(Ub-Ua))/rhob)
//                     + avCount*dragContbMean)/(avCount+1);

    if (runTime.outputTime())
    {
        
        TKEMeanProd.write();
        turbDiffMean.write();
        pressDiffMean.write();
        SGSDiffMean.write();
        viscDiffusionMean.write();
        viscDissMean.write();
        SGSDissMean.write();
        dragContbMean.write();
        // take the divergence of the time-averaged turbulent fluxes
        turbDiffusionMean = fvc::div(turbDiffMean);
        pressDiffusionMean = fvc::div(pressDiffMean);
        SGSDiffusionMean = fvc::div(SGSDiffMean);

        turbDiffusionMean.write();
        pressDiffusionMean.write();
        SGSDiffusionMean.write();

    }
}

avCount++;
