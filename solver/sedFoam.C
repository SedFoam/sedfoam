/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    SedFoam

Description
    Solver for a system of 2 phases with one phase dispersed,
    e.g. gas bubbles in a liquid or solid particles in a gas.

\*---------------------------------------------------------------------------*/
/**
 * \file sedFoam.C
 * \brief 2 phases Solver
 * \author Julien Chauchat, Cyrille Bonamy, Antoine Mathieu, Tim Nagel,
           Zhen Cheng and Tian-Jian Hsu.
 * \version 3.1
 * \date September 16, 2019
 *
 * Solver for a system of 2 phases with one phase dispersed
 *
 */
/*
* Changelog [Higuera]
* December 10, 2018 - Adapted to compile automatically with OpenFOAM 6.
*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "PhaseIncompressibleTurbulenceModel.H"

#include "symmetryFvPatchFields.H"
#include "fixedFluxPressureFvPatchScalarField.H"

#include "dragModel.H"
#include "phaseModel.H"

#include "kineticTheoryModel.H"
#include "granularRheologyModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fixedValueFvsPatchFields.H"
//#include "IOMRFZoneList.H"
//#include "IOMRFZoneList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
 //   #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"

//    #include "readGravitationalAcceleration.H"
    #include "readGravity.H"
    #include "createGradP.H"
    #include "createFields.H"
    #include "createRASTurbulence.H"
    #include "createFvOptions.H"

    #include "readPPProperties.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
 //   pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    if (SUSlocal)
    {
        Info<< "\nLocal Schmidt number activated" << endl;
    }
    else
    {
       if (max(SUS).value() == 0)
       {
           Info<< "Turbulence suspension term is neglected" << endl;
       }
       else if (max(SUS).value() > 0)
       {
           Info<< "Turbulence suspension term is included" << endl;
       }
       else
       {
           Info<< "Turbulence suspension coefficient SUS can't be negative"
               << endl;
       }
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTwoPhaseEulerFoamControls.H"
        #include "CourantNos.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "gravityRamp.H"

//      Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaEqn.H"
            #include "liftDragCoeffs.H"

//          Compute the Kinetic Theory parameters: nuEffa and lambdaUa from the
//          solution of the Granular temperature equation
            #include "callKineticTheory.H"

//          Compute the contact pressure pff and the Frictional stress nuFra
//          from a Coulomb model if using the kinetic theory
//          and from the mu(I) rheology if using the granular rheology
            #include "callFrictionStress.H"

//          Create the momentum balance equations for both phases a and b
            #include "UEqns.H"

            #include "pEqn.H"

            #include "DDtU.H"

            if (pimple.turbCorr())
            {
                #include "updateTwoPhaseRASTurbulence.H"
                turbulenceb->correct();
                if (debugInfo)
                {
                    Info << " max(nutb) = "
                         << max(turbulenceb->nut()).value() << endl;
                }
            }
        }
        if (debugInfo)
        {
            Info<< "min(Ua) = " << gMin(Ua)
                << "max(Ua) = " << gMax(Ua) << endl;
            Info<< "min(Ub) = " << gMin(Ub)
                << "max(Ub) = " << gMax(Ub) << nl << endl;
        }
        #include "OutputGradPOSC.H"
        #include "writeTau.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
