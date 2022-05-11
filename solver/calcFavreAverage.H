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
    TKEMean = 0.5 * (
                       UbPrime2MeanF.component(0)
                     + UbPrime2MeanF.component(3)
                     + UbPrime2MeanF.component(5)
                     );
    TKEMeanProd = -UbPrime2MeanF && fvc::grad(UbMeanF);
    BB = - beta*(2.0*symm(fvc::grad(Ub))
				-((2.0/3.0)*I)*tr(fvc::grad(Ub)().T()));

    turbDiffMean = (turbDiffMean 
				- avCount*0.5*(UbPrimeF*magSqr(UbPrimeF)))/(avCount+1);
    turbDiffusionMean = fvc::div(turbDiffMean);
    
    //Pressure Diffusion Term
    pressDiffusionMean = (pressDiffusionMean 
	- avCount*1./phaseb->rho()*(UbPrimeF & fvc::grad(pPrimeF)))/(avCount+1);
    SGSDiffMean = (SGSDiffMean -avCount* 
    (UbPrimeF & (turbulenceb->nut()*dev(twoSymm(fvc::grad(Ub))))))/(avCount+1);
    SGSDiffusionMean = fvc::div(SGSDiffMean);

    //Tensor of Strain Rate of Resolved Fluctuations
    SGSstrainTensor = symm(fvc::grad(UbPrimeF)); 
    
    //Viscous Dissipation of Resolved Fluctuations
    viscDissMean = (viscDissMean - avCount*2*phaseb->nu()*
    (SGSstrainTensor && SGSstrainTensor))/(avCount+1); 
    SGSDissMean = (SGSDissMean + 
    avCount*(turbulenceb->nut()*BB && SGSstrainTensor))/(avCount+1);

    dragContbMean = ((-UbPrimeF & (alpha*beta*draga->K(mag(Ub-Ua))*(Ub
					 -Ua))/rhob)
                     + avCount*dragContbMean)/(avCount+1);
//    dragContbMean = ((-UbPrimeF & (alpha*beta*(Ub-Ua))/rhob)
//                     + avCount*dragContbMean)/(avCount+1);

    if (runTime.outputTime())
    {
        TKEMeanProd.write();
        turbDiffusionMean.write();
        pressDiffusionMean.write();
        SGSDiffusionMean.write();
        viscDissMean.write();
        SGSDissMean.write();
        dragContbMean.write();
    }
}

avCount++;
