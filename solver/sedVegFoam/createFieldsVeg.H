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

Info<< "Reading vegetation properties\n" << endl;

autoPtr<phaseModel> phasec = phaseModel::New
(
    mesh,
    transportProperties,
    "c"
);

volScalarField& alphaVeg = phasec->alpha();
Info<< "alpha_b is corrected to 1 - alpha_a - alpha_c\n" << endl;
beta = 1 - alpha - alphaVeg;

volScalarField dc = phasec->d()*Unity;
volScalarField dp = phasea->d()*Unity;
// Pore space between vegetation cylinders
//volScalarField sn = min(
//      dc*(sqrt(constant::mathematical::pi/(4*(alphaVeg+alphaSmall)))-1),
//      dp*(pow(alphaMax/(alpha+alphaSmall), 1.0/3.0)-1));
volScalarField sn = dc*(
        sqrt(constant::mathematical::pi/(4*(1-beta+alphaSmall)))-1);

autoPtr<vegDragModel> dragbveg = vegDragModel::New
(
    interfacialProperties,
    phaseb,
    phasec
);

volScalarField Kveg
(
    IOobject
    (
        "Kveg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    dragbveg->K(mag(Ub))
);

volVectorField VegDragContr
(
    IOobject
    (
        "VegDragContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);
