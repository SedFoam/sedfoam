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
    e.g.  solid particles in a fluid.

    Reference:
    \verbatim
        Chauchat J., Cheng Z., Nagel T., Bonamy C., & Hsu T-J. (2017).
        SedFoam-2.0: a 3D two-phase flow numerical model for sediment transport
        Geosci. Model Dev. Discuss.
        http://dx.doi.org/10.5194/gmd-2017-101
    \endverbatim

    This file is inside the \dir solver directory

Version
    3.2

Author
    Julien Chauchat, Cyrille Bonamy, Antoine Mathieu, Rémi Chassagne,
    Tim Nagel, Zhen Cheng, Tian-Jian Hsu and Eduard Puig Montella.

Date
    June 01, 2021

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "sedIncompressibleTurbulenceModel.H"

#include "symmetryFvPatchFields.H"
#include "fixedFluxPressureFvPatchScalarField.H"

#include "dragModel.H"
#include "vegDragModel.H"
#include "phaseModel.H"
#include "ppModel.H"

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


    #include "readGravity.H"
    #include "createFields.H"
    #include "createFieldsVeg.H"
    #include "createTwoPhaseTurbulenceVeg.H"
    #include "createTurbulence.H"
    #include "createFvOptions.H"

    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "createFavreAveraging.H"
    #include "createFavreAverage_fluid.H"
    #include "createFavreAverage_particle.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Test on SUSlocal
    //
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
    // Test on granular stress model
    if (kineticTheory.on() && granularRheology.on())
    {
        Info<< "\nKinetic theory and granular rheology are set on." << endl;
        Info<< " This option is not supported!" << endl;
    }
    // stress formulation
    Switch faceMomentum
    (
        pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)
    );
    Info<< "Choice for faceMomentum : "<<faceMomentum
        << endl;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTwoPhaseEulerFoamControls.H"
        #include "CourantNos.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

//      Apply a ramp in time on the gravity acceleration
        #include "gravityRamp.H"

//      Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
//          Solve for solid phase mass conservation
            #include "alphaEqn.H"

//          Compute lift and drag coefficients
            #include "liftDragCoeffs.H"

//          Compute the granular stress: pff, nuFra, nuEffa and lambdaUa
//             from Kinetic Theory of granular flows or mu(I) rheology
            #include "callGranularStress.H"

//          Assemble the momentum balance equations for both phases a and b
//          And assemble and solve the pressure poisson equation
//             and apply the velocity correction step for both phases a and b
            if (faceMomentum)
            {
                #include "pUf/UEqns.H"
                #include "pUf/pEqn.H"
            }
            else
            {
                #include "pU/UEqns.H"
                #include "pU/pEqn.H"
            }

//          Compute the phase accelerations for added mass force
            #include "DDtU.H"

            if (pimple.turbCorr())
            {
//              Solve for turbulence models
                #include "updateTwoPhaseTurbulence.H"
                #include "updateTwoPhaseTurbulenceVeg.H"
                turbulenceb->correct();
                if (turbulencePropertiesb.get<word>("simulationType")=="LES")
                {
                    spherSigmaSGSb = turbulenceb->spherSigmaSGS();
                }
                turbulencea->correct();
                if (turbulencePropertiesa.get<word>("simulationType")=="LES")
                {
                    spherSigmaSGSa = turbulencea->spherSigmaSGS();
                }

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
//      Write output
        #include "OutputGradPOSC.H"
        #include "writeOutput.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
