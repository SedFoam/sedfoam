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

Switch favreAverage_momentum
(
    runTime.controlDict().lookupOrDefault("favreAverage_momentum", false)
);

Switch Reset_MomentumBudget
(
    runTime.controlDict().lookupOrDefault("Reset_MomentumBudget", false)
);

Info<<"Manohar Momentum Budget in create File = \t"<<favreAverage_momentum<<endl;

scalar timeaverage_momentum = 0;
scalar runtime_momentum     = 0;
double time_diff_momentum = 0;
double time_integrate_momentum = 0;

//~ volVectorField DragContrVeg_b
//~ (
    //~ IOobject
    //~ (
        //~ "DragContrVeg_b",
         //~ runTime.timeName(),
         //~ mesh,
         //~ IOobject::NO_READ,
         //~ IOobject::NO_WRITE
    //~ ),
    //~ fvc::grad(pff)
//~ );

volScalarField pMech_b
(
    IOobject
    (
        "pMech_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pMech_b",
        dimPressure,
        0
    )
);

volScalarField pA
(
    IOobject
    (
        "pA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pA",
        dimPressure,
        0
    )
);

volScalarField pFF
(
    IOobject
    (
        "pFF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pFF",
        dimPressure,
        0
    )
);

volVectorField rho_mixd
(
    IOobject
    (
        "rho_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "rho_mixd",
        dimensionSet(1, -2, -1, 0, 0, 0, 0),
        vector (0, 0, 0)
    )
);

volTensorField nonl_mixd
(
    IOobject
    (
        "nonl_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "nonl_mixd",
        dimensionSet(1, -1, -2, 0, 0, 0, 0),
        tensor::zero
    )
);

volVectorField bouy_mixd
    (
        IOobject
        (
            "bouy_mixd",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "bouy_mixd",
            dimensionSet(1, -2, -2, 0, 0, 0, 0),
            vector (0, 0, 0)
        )
    );

volVectorField SolidPressureContr_b
(
    IOobject
    (
        "SolidPressureContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(p_rbgh)
);

volVectorField divTauaContr_b
(
    IOobject
    (
        "divTauaContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField ViscStressContr_b
(
    IOobject
    (
        "ViscStressContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField FrictContr_b
(
    IOobject
    (
        "FrictContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField phiRaContr_b
(
    IOobject
    (
        "phiRaContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField phiRa2Contr_b
(
    IOobject
    (
        "phiRa2Contr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField divRcaContr_b
(
    IOobject
    (
        "divRcaContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField RcaContr_b
(
    IOobject
    (
        "RcaContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField BulkContr_b
(
    IOobject
    (
        "BulkContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField gravityContr_b
(
    IOobject
    (
        "gravityContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField DragContr_b
(
    IOobject
    (
        "DragContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField SuspContr_b
(
    IOobject
    (
        "SuspContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField ExtContr_b
(
    IOobject
    (
        "ExtContr_b",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volScalarField Errosion
        (
            IOobject
                (
                    "Errosion",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
            IOobject::NO_WRITE
                ),
            mesh,
            dimensionedScalar
                (
                    "Errosion",
                        dimensionSet(0, 1, -1, 0, 0, 0, 0),
                        0.0
                )
 );

volScalarField Deposition
        (
            IOobject
                (
                    "Deposition",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
            IOobject::NO_WRITE
                ),
            mesh,
            dimensionedScalar
                (
                    "Deposition",
                        dimensionSet(0, 1, -1, 0, 0, 0, 0),
                        0.0
                )
        );

if (Reset_particle)
{

Info<<"Manohar in Reset mode = \t"<<endl;

runtime_momentum = 0;

}

else

{

    runtime_momentum = runTime.value();// 0;// this->runTime.timeName();

IOobject alpha_aMeanHeader
            (
                "alpha_aMean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            );

if (alpha_aMeanHeader.typeHeaderOk<volScalarField>(true))
{
    volScalarField alpha_aMean2
                    (
                            IOobject
                                (
                                    "alpha_aMean",
                                    runTime.timeName(),
                                    mesh,
                                    IOobject::MUST_READ,
                                    IOobject::NO_WRITE
                                ),
                            mesh
                    );
                    alpha_aMean = alpha_aMean2;

  IOobject pMech_bHeader
    (
        "pMech_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pMech_bHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField pMech2
(
    IOobject
    (
        "pMech_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pMech_b = pMech2;
    }

IOobject pAHeader
    (
        "pA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pAHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField pA2
(
    IOobject
    (
        "pA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pA = pA2;
    }

IOobject pFFHeader
    (
        "pFF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pFFHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField pFF2
(
    IOobject
    (
        "pFF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pFF = pFF2;
    }

IOobject rho_mixdHeader
    (
        "rho_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (rho_mixdHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField rho_mixd2
(
    IOobject
    (
        "rho_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
rho_mixd = rho_mixd2;
    }

IOobject nonl_mixdHeader
    (
        "nonl_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (nonl_mixdHeader.typeHeaderOk<volTensorField>(true))
    {
volTensorField nonl_mixd2
(
    IOobject
    (
        "nonl_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
nonl_mixd = nonl_mixd2;
    }

IOobject bouy_mixdHeader
    (
        "bouy_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (bouy_mixdHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField bouy_mixd2
(
    IOobject
    (
        "bouy_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
bouy_mixd = bouy_mixd2;
    }

IOobject SolidPressureContr_bHeader
    (
        "SolidPressureContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SolidPressureContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField SolidPressureContr2
(
    IOobject
    (
        "SolidPressureContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SolidPressureContr_b = SolidPressureContr2;
    }

IOobject divTauaContr_bHeader
    (
        "divTauaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (divTauaContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField divTauaContr2
(
    IOobject
    (
        "divTauaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
divTauaContr_b = divTauaContr2;
    }

IOobject ViscStressContr_bHeader
    (
        "ViscStressContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (ViscStressContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField ViscStressContr2
(
    IOobject
    (
        "ViscStressContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
ViscStressContr_b = ViscStressContr2;
    }

IOobject FrictContr_bHeader
    (
        "FrictContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (FrictContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField FrictContr2
(
    IOobject
    (
        "FrictContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
FrictContr_b = FrictContr2;
    }

IOobject phiRaContr_bHeader
    (
        "phiRaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (phiRaContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField phiRaContr2
(
    IOobject
    (
        "phiRaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
phiRaContr_b = phiRaContr2;
    }

IOobject phiRa2Contr_bHeader
    (
        "phiRa2Contr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (phiRa2Contr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField phiRa2Contr2
(
    IOobject
    (
        "phiRa2Contr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
phiRa2Contr_b = phiRa2Contr2;
    }

IOobject divRcaContr_bHeader
    (
        "divRcaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (divRcaContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField divRcaContr2
(
    IOobject
    (
        "divRcaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
divRcaContr_b = divRcaContr2;
    }

IOobject RcaContr_bHeader
    (
        "RcaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (RcaContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField RcaContr2
(
    IOobject
    (
        "RcaContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
RcaContr_b = RcaContr2;
    }

IOobject BulkContr_bHeader
    (
        "BulkContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (BulkContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField BulkContr2
(
    IOobject
    (
        "BulkContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
BulkContr_b = BulkContr2;
    }

IOobject gravityContr_bHeader
    (
        "gravityContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (gravityContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField gravityContr2
(
    IOobject
    (
        "gravityContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
gravityContr_b = gravityContr2;
    }

IOobject DragContr_bHeader
    (
        "DragContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (DragContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField DragContr2
(
    IOobject
    (
        "DragContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
DragContr_b = DragContr2;
    }

IOobject SuspContr_bHeader
    (
        "SuspContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SuspContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField SuspContr2
(
    IOobject
    (
        "SuspContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SuspContr_b = SuspContr2;
    }

IOobject ExtContr_bHeader
    (
        "ExtContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (ExtContr_bHeader.typeHeaderOk<volVectorField>(true))
    {
volVectorField ExtContr2
(
    IOobject
    (
        "ExtContr_b",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
ExtContr_b = ExtContr2;
    }

IOobject ErrosionHeader
    (
        "Errosion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (ErrosionHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField Errosion2
(
    IOobject
    (
        "Errosion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Errosion = Errosion2;
    }

IOobject DepositionHeader
    (
        "Deposition",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (DepositionHeader.typeHeaderOk<volScalarField>(true))
    {
volScalarField Deposition2
(
    IOobject
    (
        "Deposition",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Deposition = Deposition2;
    }

/////////////////

IOobject TauaHeader
    (
        "Taua",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (TauaHeader.typeHeaderOk<volTensorField>(true))
    {
volTensorField Taua2
(
    IOobject
    (
        "Taua",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Taua = Taua2;
    }

IOobject TaubHeader
    (
        "Taub",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (TaubHeader.typeHeaderOk<volTensorField>(true))
    {
volTensorField Taub2
(
    IOobject
    (
        "Taub",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Taub = Taub2;
    }

}
}
