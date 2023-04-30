{
    if (mesh.changing())
    {
		
		
		surfaceScalarField alphaf(fvc::interpolate(alpha));
		surfaceScalarField betaf(scalar(1.0) - alphaf);
			
		
		
        volVectorField::Boundary& bafld = Ua.boundaryFieldRef();
        forAll(bafld, patchi)
        {
            if (bafld[patchi].fixesValue())
            {
                bafld[patchi].initEvaluate();
            }
        }
        volVectorField::Boundary& bbfld = Ub.boundaryFieldRef();
        forAll(bbfld, patchi)
        {
            if (bbfld[patchi].fixesValue())
            {
                bbfld[patchi].initEvaluate();
            }
        }
        
        surfaceScalarField::Boundary& phiaBfld = phia.boundaryFieldRef();
        forAll(bafld, patchi)
        {
            if (bafld[patchi].fixesValue())
            {
                bafld[patchi].evaluate();

                phiaBfld[patchi] =
                    bafld[patchi]
                  & mesh.Sf().boundaryField()[patchi];
            }
        }
        surfaceScalarField::Boundary& phibBfld = phib.boundaryFieldRef();
        forAll(bbfld, patchi)
        {
            if (bbfld[patchi].fixesValue())
            {
                bbfld[patchi].evaluate();

                phibBfld[patchi] =
                    bbfld[patchi]
                  & mesh.Sf().boundaryField()[patchi];
            }
        }
        
    }

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
        fvc::makeRelative(phia, Ua);
        adjustPhi(phia, Ua, pcorr);
        fvc::makeAbsolute(phia, Ua);
        
        fvc::makeRelative(phib, Ub);
        adjustPhi(phib, Ub, pcorr);
        fvc::makeAbsolute(phib, Ub);
    }

    mesh.setFluxRequired(pcorr.name());

    dimensionedScalar rUaAf("rUaAf", dimTime/rho.dimensions(), 1.0);
    dimensionedScalar rUbAf("rUbAf", dimTime/rho.dimensions(), 1.0);

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

        // Poisson equation for the pressure
        surfaceScalarField Dp
        (
            "(rho*(1|A(U)))",
            alphaf*rUaAf/rhoa + betaf*rUbAf/rhob
        );
        phi = alphaf*phia + betaf*phib;

        fvScalarMatrix pcorrEqn
        (        
            fvm::laplacian(Dp, pcorr) == fvc::div(phi)
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
