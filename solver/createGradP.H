/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/
Info<< "\nReading forceProperties for gradP\n" <<endl;

IOdictionary forceProperties
(
    IOobject
    (
        "forceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar tilt
(
    forceProperties.lookupOrDefault
    (
        "tilt",
        dimensionedScalar("tilt", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
    )
);

dimensionedVector gradPAMP1
(
    forceProperties.lookupOrDefault
    (
        "gradPAMP1",
        dimensionedVector("gradPAMP1",
                          dimensionSet(1, -2, -2, 0, 0, 0, 0),
                          vector(0, 0, 0) )
    )
);

dimensionedVector gradPAMP2
(
    forceProperties.lookupOrDefault
    (
        "gradPAMP2",
        dimensionedVector("gradPAMP2",
                          dimensionSet(1, -2, -2, 0, 0, 0, 0),
                          vector(0, 0, 0) )
    )
);

dimensionedVector gradPMEAN
(
    forceProperties.lookup("gradPMEAN")
);

dimensionedVector gradPOSC
(
    forceProperties.lookupOrDefault
    (
        "gradPOSC",
        dimensionedVector("gradPOSC",
                          dimensionSet(1, -2, -2, 0, 0, 0, 0),
                          vector(0, 0, 0) )
    )
);

dimensionedScalar oscpT
(
    forceProperties.lookupOrDefault
    (
        "oscpT",
        dimensionedScalar("oscpT",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          0)
    )
);

// to avoid division by 0 in gradP impose 1 if 0 is specified
dimensionedScalar oscpTMin ("oscpTMin", oscpT.dimensions(), 1e-6);
oscpT=max(oscpTMin, oscpT);

dimensionedScalar Cvm
(
    forceProperties.lookup("Cvm")
);

dimensionedScalar Cl
(
    forceProperties.lookup("Cl")
);

dimensionedScalar Ct
(
    forceProperties.lookup("Ct")
);
dimensionedScalar initTheta
(
    forceProperties.lookupOrDefault
    (
        "initTheta",
        dimensionedScalar("initTheta", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
    )
);
dimensionedScalar InitFreeze
(
    forceProperties.lookupOrDefault
    (
        "InitFreeze",
        dimensionedScalar("InitFreeze", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
    )
);
dimensionedScalar fixbeddepth
(
    forceProperties.lookupOrDefault
    (
        "fixbeddepth",
        dimensionedScalar("fixbeddepth", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0)
    )
);
dimensionedScalar ClipUa
(
    forceProperties.lookupOrDefault
    (
        "ClipUa",
        dimensionedScalar("ClipUa", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0)
    )
);

// Debuging info
Switch debugInfo
(
    forceProperties.lookupOrDefault("debugInfo", false)
);
// Sponge Layer
Switch spongeLayer
(
    forceProperties.lookupOrDefault("spongeLayer", false)
);
dimensionedScalar xSmin
(
    forceProperties.lookupOrDefault
    (
        "xSmin",
        dimensionedScalar("xSmin", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0)
    )
);
dimensionedScalar xSmax
(
    forceProperties.lookupOrDefault
    (
        "xSmax",
        dimensionedScalar("xSmax", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0)
    )
);
Info<< "Initializing with specified pressure gradient:"
<< gradPOSC.value() << endl << endl;
// open files to write gradP info
OFstream gradPOSCFile
(
    "gradPOSC.txt"
);
