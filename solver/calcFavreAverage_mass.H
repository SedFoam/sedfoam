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

Info<<"Mass runTime before the addition = \t"<<runtime_mass<<endl;


time_diff_mass = runtime_mass - T0;

//Info<<" Time diff in the favre fluid average file = \t"<<time_diff_mass<<"\t runtime_mass \t"<<runtime_mass<<"\t T0 = \t"<<T0<<endl;

if (time_diff_mass<0)
{
	time_integrate_mass = 0;
}
else
{
	time_integrate_mass = time_diff_mass;
}




timeaverage_mass = time_integrate_mass  + deltaT;


Info<<" Mass time diff mass = \t"<<time_diff_mass<<"\t runtime_mass \t"<<runtime_mass<<"\t T0 = \t"<<T0<<"\t delta T = \t"<<deltaT<<endl;
Info<<" Mass time integrate mass = \t"<<time_integrate_mass<<"\t timeaverage_mass \t"<<timeaverage_mass<<"\t deltavalue = \t"<<runTime.deltaTValue()<<endl;


alphaF = alpha - alpha_aMean;

Mean_alphaUa     = (time_integrate_mass* Mean_alphaUa  + alpha* Ua* deltaT) / (timeaverage_mass);
Mean_UaFalphaF   = (time_integrate_mass* Mean_UaFalphaF  + alphaF* UaPrimeF* deltaT) / (timeaverage_mass);
alphaFMean       = (time_integrate_mass* alphaFMean  + alphaF* deltaT) / (timeaverage_mass);
MeanUaF          = (time_integrate_mass* MeanUaF  + UaPrimeF* deltaT) / (timeaverage_mass);


if (runTime.outputTime())
{
	alphaF.write();
	alphaFMean.write();
	Mean_alphaUa.write();
	Mean_UaFalphaF.write();
	MeanUaF.write();
}

runtime_mass +=runTime.deltaTValue();
