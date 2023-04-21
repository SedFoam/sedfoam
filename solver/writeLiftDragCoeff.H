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


#include <iostream>
#include <fstream>
//using namespace std;


//string path = "ForceFiles/ForceTotFile.txt";


//std::ifstream fin(path);

//std::ofstream inFileTot;

//if(fin.fail()){
	//Info<< "creating ForceTotFile.txt "  << endl;
	
	//inFileTot.open("ForceFiles/ForceTotFile.txt", std::ios_base::app);
	//inFileTot <<  "Time" << tab << "FCont_x" << tab << "FCont_y" << tab << "FCont_z" << tab << "FCol_x" << tab << "FCol_y" << tab << "FCol_z" << tab  << "FShear_x" << tab << "FShear_y" << tab << "FShear_z"<< "\n";
	//inFileTot.close ();
//} 

label Patchid=mesh.boundaryMesh().findPatchID("cylinder");

if (runTime.writeTime()){
	
	
		string strCont1 = "ForceFiles/ForceContFile";
		string strCont2 = strCont1.append(std::to_string(runTime.value()));
		string strCont3 = strCont2.append(".txt");

		string strCol1 = "ForceFiles/ForceColFile";
		string strCol2 = strCol1.append(std::to_string(runTime.value()));
		string strCol3 = strCol2.append(".txt");

		string strShear1 = "ForceFiles/ForceShearFile";
		string strShear2 = strShear1.append(std::to_string(runTime.value()));
		string strShear3 = strShear2.append(".txt");



		string strShearFluid1 = "ForceFiles/ForceShearFluidFile";
		string strShearFluid2 = strShearFluid1.append(std::to_string(runTime.value()));
		string strShearFluid3 = strShearFluid2.append(".txt");

		std::ofstream inFileForces;
std::ifstream file(strShearFluid3);
if (!file)
{

			if (Pstream::master())
			{

		inFileForces.open(strCont3, std::ios_base::app);
		inFileForces << "x_pos" << tab << "y_pos" << tab << "z_pos" << tab << "f_con_x" << tab << "f_con_y" << tab<< "f_con_z" << tab << "vf_x" << tab << "vf_y" << tab<< "vf_z" << "\n";
		inFileForces.close();  

		inFileForces.open(strCol3, std::ios_base::app);
		inFileForces << "x_pos" << tab << "y_pos" << tab << "z_pos" << tab << "f_col_x" << tab << "f_col_y" << tab<< "f_col_z" << tab << "vf_x" << tab << "vf_y" << tab<< "vf_z" << "\n";
		inFileForces.close();  

		inFileForces.open(strShear3, std::ios_base::app);
		inFileForces << "x_pos" << tab << "y_pos" << tab << "z_pos" << tab << "f_shear_x" << tab << "f_shear_y" << tab<< "f_shear_z" << tab << "vf_x" << tab << "vf_y" << tab<< "vf_z" << "\n";
		inFileForces.close();  


		inFileForces.open(strShearFluid3, std::ios_base::app);
		inFileForces << "x_pos" << tab << "y_pos" << tab << "z_pos" << tab << "f_shearF_x" << tab << "f_shearF_y" << tab<< "f_shearF_z" << tab << "vf_x" << tab << "vf_y" << tab<< "vf_z" << "\n";
		inFileForces.close();  

	
			}
			    // file is not open
}
			
        const surfaceVectorField::Boundary& Sfb =
            mesh.Sf().boundaryField();
            
		const vectorField fN_colisional
		(
		   Sfb[Patchid]
		   *(pa.boundaryField()[Patchid])
		);

		const vectorField fN_contact
		(
		   Sfb[Patchid]
		   *(pff.boundaryField()[Patchid])
		);
		
		
		//volTensorField TauaMod = Taua;
		tensorField TauaMod = Taua.boundaryField()[Patchid].T();

		//forAll(TauaMod,cellI)
		//{
			//TauaMod[cellI].xx()=0;
			//TauaMod[cellI].yy()=0;

			////Taua.boundaryField()[Patchid][cellI][4]=0;
			//// Info << "TauaMod.boundaryField()[Patchid].xx()[cellI]: " << TauaMod[cellI]<< endl;

		//}
		

	   //TauaMod.component(tensor::XX) = 0;
	   //TauaMod.component(tensor::YY) = 0;



		const vectorField fT_shear(Sfb[Patchid] & TauaMod);


       //     Info << "Taua: " << TauaMod.boundaryField()[Patchid]<< endl;
      //      Info << "fT_shear[i]: " << fT_shear<< endl;
            
            
		const vectorField fT_shearFluid(Sfb[Patchid] & Taub.boundaryField()[Patchid]);

		const vectorField Coordinates
		(
		   mesh.C().boundaryField()[Patchid]
		);
		
        vector FContTot(0,0,0);
        vector FColiTot(0,0,0);
        vector FShearTot(0,0,0);

		
        for (int i = 0; i < fN_contact.size(); i++) {
            FContTot[0] += fN_contact[i][0];
            FContTot[1] += fN_contact[i][1];
            FContTot[2] += fN_contact[i][2];

            FColiTot[0] += fN_colisional[i][0];
            FColiTot[1] += fN_colisional[i][1];
            FColiTot[2] += fN_colisional[i][2];
            
            FShearTot[0] += fT_shear[i][0];
            FShearTot[1] += fT_shear[i][1];
            FShearTot[2] += fT_shear[i][2];
            
            //Info << "Taua: " << TauaMod.boundaryField()[Patchid][i] << endl;
            //Info << "fT_shear[i]: " << fT_shear[i]<< endl;





			// Only if the processor is the master, write to the file
			//if (Pstream::master())
			//{

				inFileForces.open(strCont3, std::ios_base::app);
				inFileForces << Coordinates[i][0] << tab << Coordinates[i][1] << tab << Coordinates[i][2] << tab << fN_contact[i][0] << tab << fN_contact[i][1] << tab<< fN_contact[i][2] << tab << Sfb[Patchid][i][0] << tab << Sfb[Patchid][i][1] << tab<< Sfb[Patchid][i][2] << "\n";
				inFileForces.close();  
				
				inFileForces.open(strCol3, std::ios_base::app);
				inFileForces << Coordinates[i][0] << tab << Coordinates[i][1] << tab << Coordinates[i][2] << tab << fN_colisional[i][0] << tab << fN_colisional[i][1] << tab<< fN_colisional[i][2] << tab << Sfb[Patchid][i][0] << tab << Sfb[Patchid][i][1] << tab<< Sfb[Patchid][i][2] << "\n";
				inFileForces.close();   
				
				inFileForces.open(strShear3, std::ios_base::app);
				inFileForces << Coordinates[i][0] << tab << Coordinates[i][1] << tab << Coordinates[i][2] << tab << fT_shear[i][0] << tab << fT_shear[i][1] << tab<< fT_shear[i][2] << tab << Sfb[Patchid][i][0] << tab << Sfb[Patchid][i][1] << tab<< Sfb[Patchid][i][2] << "\n";
				inFileForces.close();   
				
				
				inFileForces.open(strShearFluid3, std::ios_base::app);
				inFileForces << Coordinates[i][0] << tab << Coordinates[i][1] << tab << Coordinates[i][2] << tab << fT_shearFluid[i][0] << tab << fT_shearFluid[i][1] << tab<< fT_shearFluid[i][2] << tab << Sfb[Patchid][i][0] << tab << Sfb[Patchid][i][1] << tab<< Sfb[Patchid][i][2] << "\n";
				inFileForces.close(); 
			//}







		
		  }
	

	//inFileTot.open("ForceFiles/ForceTotFile.txt", std::ios_base::app);
	//inFileTot << runTime.value() << tab << FContTot[0] << tab << FContTot[1] << tab << FContTot[2] << tab << FColiTot[0]  << tab << FColiTot[1] << tab << FColiTot[2]  << tab  << FShearTot[0] << tab << FShearTot[1] << tab << FShearTot[2]<< "\n";
	//inFileTot.close();
   
	
}

