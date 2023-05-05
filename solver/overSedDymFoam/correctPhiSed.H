{
    if (mesh.changing())
    {
        volVectorField::Boundary& bfld = U.boundaryFieldRef();
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

                phiBfld[patchi] =
                    bfld[patchi]
                  & mesh.Sf().boundaryField()[patchi];
            }
        }
    }
    
    surfaceScalarField alphaf(fvc::interpolate(alpha));
    surfaceScalarField betaf(scalar(1.0) - alphaf);
	
    wordList pcorrTypes
    (
        p_rbgh.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    for (label i=0; i<p_rbgh.boundaryField().size(); i++)
    {
        if (p_rbgh.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
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
        dimensionedScalar(p_rbgh.dimensions(), Zero),
        pcorrTypes
    );

    if (pcorr.needReference())
    {
        fvc::makeRelative(phi, U);
        adjustPhi(phi, U, pcorr);
        fvc::makeAbsolute(phi, U);
        
        fvc::makeRelative(phia, Ua);
        adjustPhi(phia, Ua, pcorr);
        fvc::makeAbsolute(phia, Ua);
        
        fvc::makeRelative(phib, Ub);
        adjustPhi(phib, Ub, pcorr);
        fvc::makeAbsolute(phib, Ub);    
    }

    mesh.setFluxRequired(pcorr.name());

    dimensionedScalar rAUf("rAUf", dimTime/rho.dimensions(), 1.0);

    const cellCellStencilObject& overlap = Stencil::New(mesh);
    const labelList& cellTypes = overlap.cellTypes();
    const labelIOList& zoneIDs = overlap.zoneID();

    while (pimple.correctNonOrthogonal())
    {
	surfaceScalarField phir("phir",phia - phib);

		
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
            
			// Compute pressure gradient at cell faces
			surfaceScalarField SfGradp(pcorrEqn.flux()/rAUf);
			// Compute corrected velocity fluxes at cell faces
			//       due to pressure gradient
			surfaceScalarField phias(phia - rAUf*SfGradp);
			surfaceScalarField phibs(phib - rAUf*SfGradp);
            
			// Compute corrected relative velocity fluxes 
			surfaceScalarField phir(phias - phibs );
            
            phia = phi + betaf*(phir);
            phib = phi - alphaf*(phir);
            
        }
    }
}
