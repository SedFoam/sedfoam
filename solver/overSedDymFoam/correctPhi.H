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
if (mesh.changing())
{
    volVectorField::Boundary& bfld = Ub.boundaryFieldRef();
    forAll(bfld, patchi)
    {
        if (bfld[patchi].fixesValue())
        {
            bfld[patchi].initEvaluate();
        }
    }

    surfaceScalarField::Boundary& phiBfld = phi.boundaryFieldRef();
    forAll(bfld, patchi)
    {
        if (bfld[patchi].fixesValue())
        {
            bfld[patchi].evaluate();

            phiBfld[patchi] = bfld[patchi] & mesh.Sf().boundaryField()[patchi];
        }
    }
}
    // Initialize BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p_rbgh.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p_rbgh.boundaryField(), patchi)
    {
        if (p_rbgh.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(p.dimensions(), Zero),
        pcorrTypes
    );

{
    dimensionedScalar rAUf("rAUf", dimTime/rho.dimensions(), 1.0);

    const cellCellStencilObject& overlap = Stencil::New(mesh);
    const labelList& cellTypes = overlap.cellTypes();
    const labelIOList& zoneIDs = overlap.zoneID();

    while (pimple.correctNonOrthogonal())
    {
        label nZones = gMax(zoneIDs)+1;

        //label refCellI2 = -1;
        labelList refCells(nZones, -1);
        labelList refZones(nZones, -1);

        forAll(zoneIDs, cellI)
        {
            label zoneId = zoneIDs[cellI];
            if
            (
                refCells[zoneId] == -1
             && cellTypes[cellI] == cellCellStencil::CALCULATED
             && refZones[zoneId] == -1
            )
            {
                refCells[zoneId] = cellI;
                refZones[zoneId] = zoneId;
            }
        }

        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi)
        );


        // Only set reference for cells that are CALCULATED
        {
            DynamicList<label> validCells(refCells.size());
            forAll(refCells, zoneId)
            {
                if (refCells[zoneId] != -1)
                {
                    validCells.append(refCells[zoneId]);
                }
            }

            pcorrEqn.setReferences
            (
                validCells,
                scalar(0),
                true
            );
        }

        const dictionary& d = mesh.solver
        (
            pcorr.select
            (
                pimple.finalInnerIter()
            )
        );
        //Bypass virtual layer
        mesh.fvMesh::solve(pcorrEqn, d);

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }
}
